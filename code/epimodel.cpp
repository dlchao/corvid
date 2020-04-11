/* class EpiModel implementation
 *
 * Implements a stochastic coronavirus epidemic simulator.
 */

#include <cstdlib>
#include <cstring>
#include <climits>
#include <assert.h>
extern "C" {
  #include "dSFMT.h"   // for random number generation
  #include "bnldev.h"  // for binomial random numbers
}
#include "params.h"
#include "epimodel.h"
#include "epimodelparameters.h"
#include <malloc.h> 

using namespace std;

const int nVersionMajor = 0;
const int nVersionMinor = 5;

#define get_rand_double dsfmt_genrand_close_open(&dsfmt)
#define get_rand_uint32 dsfmt_genrand_uint32(&dsfmt)

// vaccine efficacy over time when individuals need a boost for a one-dose vaccine
const double defaultvacceff[VACCEFFLENGTH+1] = {0,0.001,0.004,0.011,0.023,0.043,0.07,0.106,0.153,0.211,0.28,0.363,0.46,0.572,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.702,0.71,0.73,0.766,0.82,0.897,1};
const int defaultboostday = 21;

EpiModel::EpiModel(EpiModelParameters &params) {
  cout << "Corvid version " << nVersionMajor << "." << nVersionMinor << endl;
  nLogFileInterval = params.getLogFileInterval();
  logfile = sumfile = individualsfile = NULL;
  szBaseName = params.getBaseName();
  szLabel = params.getLabel();
  seeddisp = params.getRandomSeed();
  bTractFile = true;
  bIndividualsFile = false;
  tractToFIPStract = tractToFIPScounty = tractToFIPSstate = NULL;
  nNumPerson=nNumFamilies=nNumCommunities=nFirstTract=nLastTract=0;
  for (int i=0; i<NUMVACCINES; i++)
    nNumVaccineDosesUsed[i] = 0;
  nNumAntiviralsUsed=0;
  nNumTAPDone=0;
  nNumWantAV = nNumWantVaccine = 0;
  nTriggerTime=INT_MAX;
  bTrigger=false;
  nTimer = 0;

  // epidemic parameters
  fResponseThreshold=params.getResponseThreshold();
  bSeedDaily = params.getSeedDaily();
  nSeedAirports = params.getSeedAirports();
  bTravel = params.getTravel();
  nSeedInfectedTractFIPS = params.getSeedInfectedTractFIPS();
  nSeedInfectedCountyFIPS = params.getSeedInfectedCountyFIPS();
  nSeedInfectedStateFIPS = params.getSeedInfectedStateFIPS();
  nSeedInfectedNumber = params.getSeedInfectedNumber();
  memcpy(nSchoolOpeningDays, params.getSchoolOpeningDays(), 56*sizeof(int));

  beta = params.getBeta();
  R0 = params.getR0();
  memcpy(seasonality, params.getSeasonality(), MAXRUNLENGTH*sizeof(double));

  nRunLength = params.getRunLength();
  fPreexistingImmunityProtection = params.getPreexistingImmunityProtection();
  memcpy(fPreexistingImmunityFraction, params.getPreexistingImmunityByAge(), TAG*sizeof(double));
  memcpy(fBaselineVESByAge, params.getBaselineVESByAge(), TAG*sizeof(double));

  // vaccination/AVs
  ePrevaccinationStrategy=params.getPrevaccinationStrategy();
  eVaccinationStrategy=params.getVaccinationStrategy();
  memcpy(nVaccinePriorities, params.getVaccinePriorities(), PRIORITY_LAST*sizeof(unsigned char));
  nHighestPriority = 0;
  for (int i=0; i<PRIORITY_LAST; i++)
    nHighestPriority = max(nHighestPriority, nVaccinePriorities[i]);

  memcpy(nVaccinePriorities2, params.getVaccinePriorities2(), PRIORITY_LAST*sizeof(unsigned char));
  nHighestPriority2 = 0;
  for (int i=0; i<PRIORITY_LAST; i++)
    nHighestPriority2 = max(nHighestPriority2, nVaccinePriorities2[i]);
  nPriorityChangeTime=params.getPriorityChangeTime();

  fVaccinationFraction=params.getVaccinationFraction();
  vaccinedatastruct *vdata = params.getVaccineData();
  for (int i=0; i<NUMVACCINES; i++) {
    memcpy(VaccineData+i, vdata+i, sizeof(vaccinedatastruct));
    memcpy(vaccineproductionschedule+i, params.getVaccineProductionSchedule(i), MAXRUNLENGTH*sizeof(unsigned int));
    if (VaccineData[i].nNumDoses<=0)
      break;
    nNumVaccineTypes = i+1;
  }
  memcpy(fVaccineEfficacyByAge, params.getVaccineEfficacyByAge(), TAG*sizeof(double));
  memcpy(bVaccineBoostByAge, params.getVaccineBoostByAge(), (TAG+3)*sizeof(bool));
  memcpy(nVaccineInitialSupply, params.getVaccineInitialSupply(), NUMVACCINES*sizeof(unsigned int));
  memcpy(nVaccineSupply, nVaccineInitialSupply, NUMVACCINES*sizeof(unsigned int));
  nAVTotalLimit = params.getAVTotalLimit(); 
  nVaccineDailyLimit = params.getVaccineDailyLimit(); 
  nAVDailyLimit = params.getAVDailyLimit(); 
  nTriggerDelay=params.getTriggerDelay();
  if (nTriggerDelay<0) {
    nTriggerTime=1;  // so that the response will start on the first evening of the simulation
    bTrigger=true;
  }
  int temp=params.getTriggerDay();
  if (temp>=0) {
    nTriggerTime = 2*temp+1; // force the response to occur on specified day. nTriggerTime needs to be an odd number because interventions are triggered during the day.
    bTrigger=true;
  }
  nAscertainmentDelay = params.getAscertainmentDelay();
  fSymptomaticAscertainment=params.getAscertainmentFraction();
  fContactAscertainment = params.getContactAscertainmentFraction();
  eAVPolicy=params.getAVPolicy();
  fAdultEssentialFraction = params.getAdultEssentialFraction();
  memcpy(fPregnantFraction, params.getPregnantFraction(), TAG*sizeof(double));
  memcpy(fHighRiskFraction, params.getHighRiskFraction(), TAG*sizeof(double));

  AVEs = params.getAVEs();
  AVEi = params.getAVEi();
  AVEp = params.getAVEp();
  
  // NPIs
  nSchoolClosurePolicy=params.getSchoolClosurePolicy();
  nSchoolClosureDays=params.getSchoolClosureDays();
  fVoluntaryIsolationCompliance=params.getVoluntaryIsolationCompliance();
  fAscertainedIsolationCompliance=params.getAscertainedIsolationCompliance();
  fQuarantineCompliance=params.getQuarantineCompliance();
  fLiberalLeaveCompliance=params.getLiberalLeaveCompliance();
  fWorkFromHomeCompliance=params.getWorkFromHomeCompliance();
  nLiberalLeaveDuration=params.getLiberalLeaveDuration();
  nWorkFromHomeDuration=params.getWorkFromHomeDuration();
  nQuarantineLength=params.getQuarantineLength();
  fCommunityContactReduction=params.getCommunityContactReduction();
  nCommunityContactReductionDuration=params.getCommunityContactReductionDuration();

  // make cumulative distribution for withdraw probabilities
  // convert from double to unsigned int for efficiency
  for (int j=0; j<3; j++) {
    double sum=0.0;
    withdrawcdf[j][0] = withdrawprob[j][0];
    for (int i=1; i<WITHDRAWDAYS; i++) {
      sum += withdrawprob[j][i-1];
      withdrawcdf[j][i] = withdrawcdf[j][i-1] + ((1.0-sum)*withdrawprob[j][i]);
      assert(withdrawcdf[j][i]<=1.0);
    }
  }

  // initialize RNG
  dsfmt_init_gen_rand(&dsfmt, seeddisp);

  double vmin = basevload[0][0];
  double vmax = basevload[0][0];
  for (int i = 0; i < VLOADNSUB; i++)
    for (int j = 0; j < VLOADNDAY; j++) {
      if (basevload[i][j] < vmin)
	vmin = basevload[i][j];
      if (basevload[i][j] > vmax)
	vmax = basevload[i][j];
    }

  // assuming linear relationship between 
  // viral load and infectiousness
  double scale = vmax - vmin;
  scale *= fRelativeSymptomaticInfectiousness; // because we don't want the product of vload and fRelativeSymptomaticInfectiousness to go above 1.0
  for (int i = 0; i < VLOADNSUB; i++)
    for (int j = 0; j < VLOADNDAY; j++)
      vload[i][j] = beta * (basevload[i][j] - vmin) / scale;

  {
    cout << "Parameter set: " << szLabel << "\n";
    cout << "  " << szBaseName << " population and workflow data" << endl;
    if (R0>0.0) {
      cout << "  " << R0 << " read in for R0" << endl;
      cout << "   interpolated beta is " << beta << endl;
    } else {
      cout << "  " << beta << " read in for beta" << endl;
    }
  }

  initPopulation();

  ofstream *tractfile = params.getTractFile();
  if (tractfile) {
    ostream &out = *tractfile;
    out << "TractID,FIPSstate,FIPScounty,FIPStract,pop0-4,pop5-18,pop19-29,pop30-64,pop65+,workers" << endl;

    for (vector< Tract >::iterator it = tractvec.begin();
	 it != tractvec.end();
	 it++) {
      Tract &t = *it;
      int ntot[TAG];	// number of people
      memset(ntot, 0, sizeof(int)*TAG);
      int totalpop = 0;
      int totalworkers = 0;
      for (unsigned int i=t.nFirstCommunity; i<t.nLastCommunity; i++) {
	totalpop += commvec[i].nNumResidents;
	totalworkers += commvec[i].nNumWorkers;
	for (int j=0; j<TAG; j++) {
	  ntot[j] += commvec[i].nNumAge[j];
	}
      }
      out << t.id << "," << t.fips_state << "," << t.fips_county << "," << t.fips_tract;
      for (int j=0; j<TAG; j++) {
	out << "," << ntot[j];
      }
      out << "," << totalworkers << endl;
    }

    (*tractfile).close();
  }
    
  {
    logfile = params.getLogFile();
    if ((*logfile) && (*logfile).good()) {
      ostream &out = *logfile;
      out << "time,TractID,sym0-4,sym5-18,sym19-29,sym30-64,sym65+,cumsym0-4,cumsym5-18,cumsym19-29,cumsym30-64,cumsym65+,inf0-4,inf5-18,inf19-29,inf30-64,inf65+,cuminf0-4,cuminf5-18,cuminf19-29,cuminf30-64,cuminf65+" << endl;
    }
    individualsfile = params.getIndividualsFile();
    sumfile = params.getSummaryFile();
  }
  bIndividualsFile = params.getHasIndividualsFile();

  vector< Person >::iterator pend=pvec.end();
  for (vector< Person >::iterator it = pvec.begin(); 
       it != pend;
       it++) {
    Person &p = *it;
    if (fPreexistingImmunityFraction[p.age]>0.0 && get_rand_double<fPreexistingImmunityFraction[p.age])
      p.fBaselineVES = fPreexistingImmunityProtection;
    else
      p.fBaselineVES = fBaselineVESByAge[p.age];
    if (ePrevaccinationStrategy==PRIMEBOOSTRANDOM || // everyone is eligible when PRIMEBOOSTRANDOM
	(fVaccinationFraction>0.0 && get_rand_double<fVaccinationFraction)) {
      unsigned char priority = 100;
      if (nVaccinePriorities[PRIORITY_PREGNANT]>0 && isPregnant(p))
	priority = min(nVaccinePriorities[PRIORITY_PREGNANT], priority);
      if (nVaccinePriorities[PRIORITY_ESSENTIAL]>0 && isEssential(p))
	priority = min(nVaccinePriorities[PRIORITY_ESSENTIAL], priority);
      if (nVaccinePriorities[PRIORITY_HR0+p.age]>0 && isHighRisk(p))
	priority = min(nVaccinePriorities[PRIORITY_HR0+p.age], priority);
      if (nVaccinePriorities[PRIORITY_0+p.age]>0)
	priority = min(nVaccinePriorities[PRIORITY_0+p.age], priority);
      if (nVaccinePriorities[PRIORITY_INFANTFAMILY]>0 && !isInfant(p)) {
	unsigned int upper = (p.id+7>pvec.size()?pvec.size():p.id+7);
	for (unsigned int i=(p.id>7?p.id-7:0); i<upper; i++) {
	  if (pvec[i].family==p.family &&
	      isInfant(pvec[i])) {
	    priority = min(nVaccinePriorities[PRIORITY_INFANTFAMILY], priority);
	    i=upper;
	  }
	}
      }
      if (priority>nHighestPriority)
	p.nVaccinePriority = 0;
      else
	p.nVaccinePriority = priority;
      for (int vacnum=0; vacnum<nNumVaccineTypes; vacnum++)
	p.bVaccineEligible[vacnum] = isEligible(p, vacnum);
      if (bVaccineBoostByAge[p.age])
	setNeedsBoost(p);
      if (!needsBoost(p) && p.age==1 && p.nWorkplace>0) { // schoolchild
	if ((p.nWorkplace==1 && bVaccineBoostByAge[TAG+2]) || // high school
	    (p.nWorkplace==2 && bVaccineBoostByAge[TAG+1]) || // middle school
	    ((p.nWorkplace==3 || p.nWorkplace==4) && bVaccineBoostByAge[TAG+0]))  // elementary school
	  setNeedsBoost(p);
      }
    }
  }
}

/*
 * Read in census tracts information and create communities
 */
void EpiModel::read_tracts(void) {
  ostringstream oss;
  oss.str(szBaseName+"-tracts.dat");
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << oss.str() << " not found." << endl;
    exit(-1);
  }
  istream_iterator< Tract > iscit(iss), eos;
  copy(iscit, eos, inserter(tractvec, tractvec.begin()));
  nNumTractsTotal = tractvec.size();
  tractToFIPStract = new unsigned int[nNumTractsTotal];
  tractToFIPScounty = new unsigned int[nNumTractsTotal];
  tractToFIPSstate = new unsigned int[nNumTractsTotal];
  for (unsigned int i=0; i<tractvec.size(); i++) {
    tractToFIPStract[i] = tractvec[i].fips_tract;
    tractToFIPScounty[i] = tractvec[i].fips_county;
    tractToFIPSstate[i] = tractvec[i].fips_state;
  }
  unsigned int nEstPop = 0; // estimated total population
  for (vector< Tract >::iterator it = tractvec.begin();
       it != tractvec.end();
       it++)
    nEstPop += (*it).censuspopulation;
  pvec.reserve( (unsigned int)(nEstPop*1.0005) ); // set aside estimated memory for population (with a little padding)
  commvec.reserve( (unsigned int)(nEstPop/(double)(Community::TARGETCOMMUNITYSIZE))+100 ); // set aside estimated memory for communities

  // allocate communities
  nLastTract = nFirstTract;
  for (vector< Tract >::iterator it = tractvec.begin();
       it != tractvec.end();
       it++) {
    Tract &t = *it;
    t.id=nLastTract++;
    t.status=0; // no interventions in place
    t.nFirstCommunity = nNumCommunities;
    t.nNumResidents=0;
    t.nSchoolClosureTimer=0;
    t.nLiberalLeaveTimer=0;
    t.nWorkFromHomeTimer=0;
    t.nCommunityContactReductionTimer=0;
    for (int i=0; i<9; i++)
      t.bSchoolClosed[i] = false;
    int ncom = int(t.censuspopulation/(double)(Community::TARGETCOMMUNITYSIZE)+0.5);
    // small tracts are set to have no residents
    if (t.censuspopulation<500)
      t.censuspopulation=0;
    if (ncom==0)
      ncom=1; // we might need this for daytime workers
    while (ncom > 0) {
      Community c;
      c.id=nNumCommunities++;
      c.nNumWorkers=0;
      c.nNumNonWorkers=0;
      c.nNumWorkersLeaving=0;
      c.nTractID=t.id;
      create_families(c, (t.censuspopulation-t.nNumResidents)/ncom);
      t.nNumResidents+=c.nNumResidents;
      commvec.push_back(c);
      ncom--;
    }
    t.nLastCommunity = nNumCommunities;
  }
  cout << tractvec.size() << " tracts, " << commvec.size() << " communities" << endl;
  return;
}

// Read one tract-to-tract workflow file (called by read_workflow)
bool EpiModel::read_workflow_file(string s, unsigned int *flow, map<int,int> *statetracts) {
  ifstream iss(s.c_str());
  if (!iss)
    return false;
  int nNumTracts = tractvec.size();  // number of tracts on this node
  while (iss) {
    unsigned int fromstate, fromcounty, fromtract, tostate, tocounty, totract, workers;
    if (iss >> fromstate >> fromcounty >> fromtract >> tostate >> tocounty >> totract >> workers) {
      int from=-1, to=-1;
      map<int,int>::iterator it = statetracts[fromstate].find(fromcounty*1000000+fromtract);
      if( it != statetracts[fromstate].end() ) {
	from=it->second;
      }
      it = statetracts[tostate].find(tocounty*1000000+totract);
      if ( it != statetracts[tostate].end() )
	to=it->second;
      if (from>=0 && to>=0) {
	assert(from*nNumTractsTotal + to<nNumTracts*nNumTractsTotal);
        flow[from*nNumTractsTotal + to] = workers;
      }
    }
  }
  iss.close();
  return true;
}

// Read tract-to-tract workflow information
void EpiModel::read_workflow(void)
{
  // create map of FIPS
  map<int,int> statetracts[57]; // one tract map for each state (so we don't get integer overflows)
  for (unsigned int i=0; i<nNumTractsTotal; i++) {
    // a fast way to convert a FIPS code to an index in the tract array
    statetracts[tractToFIPSstate[i]].insert(make_pair(tractToFIPScounty[i]*1000000+tractToFIPStract[i],i));
    assert(tractToFIPStract[i]<1000000);
  }
  int nNumTracts = tractvec.size();  // number of tracts on this node

  // Read employment data
  ostringstream oss;
  oss.str(szBaseName+"-employment.dat");
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << oss.str() << " not found." << endl;
    exit(-1);
  }

  // Parse employment data
  while (iss) {
    unsigned int state, county, tract, employed, workforce;
    if (iss >> state >> county >> tract >> employed >> workforce) {
      int from=-1;
      map<int,int>::iterator it = statetracts[state].find(county*1000000+tract);
      if( it != statetracts[state].end() ) {
	from=it->second;
	if (from>=0 && from<nNumTracts) {
	  tractvec[from].employed = employed;
	  tractvec[from].workforce = workforce;
	  tractvec[from].fEmploymentProb = employed/(double)workforce;
	}
      }
    }
  }
  iss.close();

  // Read workflow data
  unsigned int *flow = new unsigned int[nNumTracts*nNumTractsTotal]; // Workerflow matrix
  if (!flow) {
    cerr << "ERROR: Could not allocate workerflow matrix. (" << (nNumTracts*nNumTractsTotal) << " ints)" << endl;
    exit(-1);
  }
  memset(flow, 0, nNumTracts*nNumTractsTotal*sizeof(unsigned int));
  oss.str(szBaseName+"-wf.dat");
  if (!read_workflow_file(oss.str(), flow, statetracts)) {
    // look at state-based wf files here
    bool bSuccess = false;
    for (unsigned int i=1; i<=56; i++) {
      bool bHasState = false;
      for (vector< Tract >::iterator it = tractvec.begin();
	   it != tractvec.end();
	   it++) { // check to see if this node has this state
	Tract &t = *it;
	if (t.fips_state==i) {
	  bHasState = true;
	  break;
	}
      }
      if (bHasState) {
	char s[3]; 
	sprintf(s, "%02d", i);
	oss.str(szBaseName+"-wf-"+s+".dat");
	if (read_workflow_file(oss.str(), flow, statetracts))
	  bSuccess = true;
      }
    }
    if (!bSuccess) {
      cerr << "ERROR: Could not find any workerflow files." << endl;
      exit(-1);
    }
  }

  // Convert to workflow cumulative distribution to enable random selection
  for (int i=0; i<nNumTracts; i++)
    for (unsigned int j=1; j<nNumTractsTotal; j++)
      flow[i*nNumTractsTotal + j] += flow[i*nNumTractsTotal + j-1];

  // Assign workplaces
  vector< Community >::iterator cend = commvec.end();
  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    unsigned int fromtract = comm.nTractID-nFirstTract;
    for (unsigned int pid=comm.nFirstPerson;
	 pid<comm.nLastPerson;
	 pid++) {
      Person &p = pvec[pid];
      if (isWorkingAge(p)) { // is adult?
	if (get_rand_double < tractvec[comm.nTractID-nFirstTract].fEmploymentProb) { // is employed?
	  if (fAdultEssentialFraction>0.0 && get_rand_double < fAdultEssentialFraction)
	    setEssential(p);

	  p.nWorkplace=-1; // assign placeholder workplace
	  p.nWorkNeighborhood=get_rand_uint32 % 4; // work neighborhood
	  p.nDayNeighborhood=p.nWorkNeighborhood; // daytime neighborhood
	  if (flow[fromtract*nNumTractsTotal + nNumTractsTotal-1]==0) {
	    // no workerflow data for this tract, so work in home tract
	    p.nDayTract = comm.nTractID;
	    if (tractvec[fromtract].nLastCommunity-tractvec[fromtract].nFirstCommunity>1 && (get_rand_double>=0.25)) {
	      p.nDayComm = tractvec[fromtract].nFirstCommunity + (get_rand_uint32 % (tractvec[fromtract].nLastCommunity-tractvec[fromtract].nFirstCommunity));
	      if (p.nDayComm!=p.nHomeComm)
		commvec[p.nDayComm].workers.push_back(pid); // add to other community worker list
	    } else // the probability of working in the home comm is high
	      p.nDayComm = p.nHomeComm;
	    commvec[p.nDayComm].nNumWorkers++;
	  } else {
	    // Choose a random destination unit
	    unsigned int irnd = get_rand_uint32 % flow[fromtract*nNumTractsTotal + nNumTractsTotal-1];
	    p.nDayTract=0;
	    while (irnd > flow[fromtract*nNumTractsTotal + p.nDayTract])
	      p.nDayTract++;
	    unsigned int totract = p.nDayTract;
	      if (tractvec[totract].nLastCommunity-tractvec[totract].nFirstCommunity>1) {
		if (fromtract==totract && (get_rand_double<0.25)) // the probability of working in the home comm is high
		  p.nDayComm=p.nHomeComm;
		else
		  p.nDayComm=tractvec[totract].nFirstCommunity + (get_rand_uint32 % (tractvec[totract].nLastCommunity-tractvec[totract].nFirstCommunity));
	      } else
		p.nDayComm=tractvec[totract].nFirstCommunity;
	      commvec[p.nDayComm].nNumWorkers++;
	      if (p.nDayComm!=p.nHomeComm)
		commvec[p.nDayComm].workers.push_back(pid); // add to other community worker list
	  }
	} else {
	  p.nWorkplace=0;   // unemployed
	  comm.nNumNonWorkers++;
	};
	assert(p.nDayTract>=0);
      } else {
	comm.nNumNonWorkers++;
      }
    }
  }
  delete flow;

  // compute number of work groups per community
  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    comm.nNumWorkGroups = (int)(comm.nNumWorkers/(float)Community::WORKGROUPSIZE + 0.5);
    if (comm.nNumWorkGroups==0 && comm.nNumWorkers>0)
      comm.nNumWorkGroups=1;
  }
  
  // assign work groups
  vector< Person >::iterator pend = pvec.end();
  for (vector< Person >::iterator it = pvec.begin(); 
       it != pend;
       it++) {
    Person &p = *it;
    if (isWorkingAge(p) && p.nWorkplace<0) {  // is adult and employed?
      p.nWorkplace = 1 + (get_rand_uint32 % commvec[p.nDayComm].nNumWorkGroups);
      if (p.nDayComm!=p.nHomeComm)
	commvec[p.nHomeComm].nNumWorkersLeaving++;
    }
  }
}

void EpiModel::create_person(int nAgeGroup, int nFamilySize, int nFamily, int nHouseholdCluster, int nNeighborhood, int nSchoolgroup, Community& comm) {
  Person p;
  p.id    = nNumPerson++;
  p.nHomeComm = p.nDayComm= comm.id; // assume home community = work community
  p.nDayTract = comm.nTractID;
  p.householdcluster = nHouseholdCluster;
  p.nWorkplace = nSchoolgroup;
  p.age   = nAgeGroup;
  p.family= nFamily;
  p.nFamilySize = nFamilySize;
  p.nHomeNeighborhood=p.nDayNeighborhood=p.nWorkNeighborhood=nNeighborhood;
  p.status = p.iday = p.ibits = p.vbits = p.nTravelTimer = p.sourcetype = 0;
  setSusceptible(p);
  p.sourceid = p.id;  // can't remember why I set the default sourceid = id
  p.nInfectedTime = -1;
  p.nIncubationDays = -1;
  p.nVaccinePriority=0;
  p.nWhichVload = get_rand_uint32%VLOADNSUB;
  p.nVaccineRestrictionBits = 0;
  memset(p.bVaccineEligible, 0, NUMVACCINES*sizeof(bool));
  if (p.age==0 && get_rand_double<fAG0InfantFraction)
    setInfant(p);
  if (fHighRiskFraction[p.age]>0.0 && (fHighRiskFraction[p.age]==1.0 || get_rand_double < fHighRiskFraction[p.age]))
    setHighRisk(p);
  if (fPregnantFraction[p.age]>0.0 && get_rand_double<fPregnantFraction[p.age])
    setPregnant(p);

  p.bWantAV = p.bWantVac = false;
  pvec.push_back(p);
}

/*
 * Fill a community with families
 * Note: work groups are not assigned here, but schools are
 * This routine was basically copied from EpiCast
 */
void EpiModel::create_families(Community& comm, int nTargetSize) {
  int  playgroup[4], // index number of the playgroup
    nplay[4],        // number of children in the playgroup
    nPlaygroupIndex = 13; // index of the next playgroup to create
  comm.nNumResidents=0;
  memset(comm.ninf, 0, TAG*sizeof(int));
  memset(comm.nsym, 0, TAG*sizeof(int));
  memset(comm.nEverInfected, 0, TAG*sizeof(int));
  memset(comm.nEverSymptomatic, 0, TAG*sizeof(int));
  memset(comm.nEverAscertained, 0, TAG*sizeof(int));
  memset(comm.nRecentlyAscertained, 0, TAG*sizeof(int));
  memset(comm.nNumAge, 0, TAG*sizeof(int));
  comm.nFirstPerson=nNumPerson;
  playgroup[0] = 9;
  playgroup[1] = 10;
  playgroup[2] = 11;
  playgroup[3] = 12;
  nplay[0] = nplay[1] = nplay[2] = nplay[3] = 0;
  while (comm.nNumResidents < nTargetSize) {
    // Create a new family and randomly place it
    // This uses the US Census's PUMS 1% data for the continental US+DC
    // (http://www2.census.gov/census_2000/datasets/PUMS/OnePercent/).
    // The distribution of household sizes was based on the PUMS data
    // for households of 1-7 people.
    // The probablities for the actual age structure of households
    // were also used  (e.g., the number of 2 elderly people living 
    // together, or one pre-schooler and two young adults, etc).
    // For households of size 1-2, all combinations of ages were included.
    // For households of size 3-7, only the household compositions that 
    // covered at least 1% of households of that size were used.
    // So if there was an unusual structure (e.g., one elderly and several
    // preschoolers), it was not included.
    int agegroups[TAG] = {0,0,0,0,0};
    double r = get_rand_double;
    if (r<0.2606) { // single-person household
      if (r<0.0955)
	agegroups[4] = 1;
      else if (r<0.2309)
	agegroups[3] = 1;
      else if (r<0.2601)
	agegroups[2] = 1;
      else
	agegroups[1] = 1;
    } else if (r<0.5892) { // two-person household
      if (r<0.3254329) {    // two elderly
	agegroups[4] = 2;
      } else if (r<0.3623) {  // elderly+old
	agegroups[4] = 1; agegroups[3] = 1;
      } else if (r<0.5037) {  // two old
	agegroups[3] = 2;
      } else if (r<0.5051) {  // elderly+young
	agegroups[4] = 1; agegroups[2] = 1;
      } else if (r<0.5295) {  // old+young
	agegroups[3] = 1; agegroups[2] = 1;
      } else if (r<0.5594) {  // two young
	agegroups[2] = 2;
      } else if (r<0.5600) {  // elderly+school
	agegroups[4] = 1; agegroups[1] = 1;
      } else if (r<0.5797) {  // old+school
	agegroups[3] = 1; agegroups[1] = 1;
      } else if (r<0.5831) {  // young+school
	agegroups[2] = 1; agegroups[1] = 1;
      } else if (r<0.5833) {  // two school
	agegroups[1] = 2;
      } else if (r<0.5853) {  // old+preschool
	agegroups[3] = 1; agegroups[0] = 1;
      } else if (r<0.5891) {  // young+preschool
	agegroups[2] = 1; agegroups[0] = 1;
      } else {                // school+preschool
	agegroups[1] = 1; agegroups[0] = 1;
      }
    } else if (r<0.7553) { // three-person household
      if (r<0.59592) {
	agegroups[3]=1;agegroups[4]=2;
      } else if (r<0.6033) {
	agegroups[3]=2;agegroups[4]=1;
      } else if (r<0.61044) {
	agegroups[3]=3;
      } else if (r<0.61275) {
	agegroups[2]=1;agegroups[3]=1;agegroups[4]=1;
      } else if (r<0.63532) {
	agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.63874) {
	agegroups[2]=2;agegroups[3]=1;
      } else if (r<0.64236) {
	agegroups[2]=3;
      } else if (r<0.64464) {
	agegroups[1]=1;agegroups[3]=1;agegroups[4]=1;
      } else if (r<0.69282) {
	agegroups[1]=1;agegroups[3]=2;
      } else if (r<0.70003) {
	agegroups[1]=1;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.70258) {
	agegroups[1]=1;agegroups[2]=2;
      } else if (r<0.71606) {
	agegroups[1]=2;agegroups[3]=1;
      } else if (r<0.73057) {
	agegroups[0]=1;agegroups[3]=2;
      } else if (r<0.73841) {
	agegroups[0]=1;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.75046) {
	agegroups[0]=1;agegroups[2]=2;
      } else if (r<0.75293) {
	agegroups[0]=1;agegroups[1]=1;agegroups[3]=1;
      } else if (r<0.75532) {
	agegroups[0]=1;agegroups[1]=1;agegroups[2]=1;
      }
    } else if (r<0.8988042) { // four-person household
      if (r<0.76155) {
	agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.76334) {
	agegroups[1]=1;agegroups[3]=2;agegroups[4]=1;
      } else if (r<0.76525) {
	agegroups[1]=1;agegroups[3]=3;
      } else if (r<0.77826) {
	agegroups[1]=1;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.84032) {
	agegroups[1]=2;agegroups[3]=2;
      } else if (r<0.84477) {
	agegroups[1]=2;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.84981) {
	agegroups[1]=3;agegroups[3]=1;
      } else if (r<0.85167) {
	agegroups[0]=1;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.87035) {
	agegroups[0]=1;agegroups[1]=1;agegroups[3]=2;
      } else if (r<0.8757) {
	agegroups[0]=1;agegroups[1]=1;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.8801) {
	agegroups[0]=1;agegroups[1]=1;agegroups[2]=2;
      } else if (r<0.88185) {
	agegroups[0]=1;agegroups[1]=2;agegroups[3]=1;
      } else if (r<0.89) {
	agegroups[0]=2;agegroups[3]=2;
      } else if (r<0.89344) {
	agegroups[0]=2;agegroups[2]=1;agegroups[3]=1;
      } else {
	agegroups[0]=2;agegroups[2]=2;
      }
    } else if (r<0.9661342) { // five-person household
      if (r<0.89985) {
	agegroups[2]=3;agegroups[3]=2;
      } else if (r<0.90282) {
	agegroups[1]=1;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.90443) {
	agegroups[1]=2;agegroups[3]=2;agegroups[4]=1;
      } else if (r<0.90618) {
	agegroups[1]=2;agegroups[3]=3;
      } else if (r<0.91267) {
	agegroups[1]=2;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.93739) {
	agegroups[1]=3;agegroups[3]=2;
      } else if (r<0.93918) {
	agegroups[1]=3;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.94066) {
	agegroups[1]=4;agegroups[3]=1;
      } else if (r<0.94171) {
	agegroups[0]=1;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.94352) {
	agegroups[0]=1;agegroups[1]=1;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.95513) {
	agegroups[0]=1;agegroups[1]=2;agegroups[3]=2;
      } else if (r<0.95772) {
	agegroups[0]=1;agegroups[1]=2;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.9589) {
	agegroups[0]=1;agegroups[1]=2;agegroups[2]=2;
      } else if (r<0.96292) {
	agegroups[0]=2;agegroups[1]=1;agegroups[3]=2;
      } else if (r<0.96458) {
	agegroups[0]=2;agegroups[1]=1;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.96613) {
	agegroups[0]=2;agegroups[1]=1;agegroups[2]=2;
      }
    } else if (r<0.9913266) { // six-person household
      if (r<0.96663) {
	agegroups[1]=1;agegroups[2]=3;agegroups[3]=2;
      } else if (r<0.96698) {
	agegroups[1]=2;agegroups[3]=4;
      } else if (r<0.96745) {
	agegroups[1]=2;agegroups[2]=1;agegroups[3]=3;
      } else if (r<0.96874) {
	agegroups[1]=2;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.96937) {
	agegroups[1]=3;agegroups[3]=2;agegroups[4]=1;
      } else if (r<0.97021) {
	agegroups[1]=3;agegroups[3]=3;
      } else if (r<0.97248) {
	agegroups[1]=3;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.97889) {
	agegroups[1]=4;agegroups[3]=2;
      } else if (r<0.97943) {
	agegroups[1]=4;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.97981) {
	agegroups[1]=5;agegroups[3]=1;
      } else if (r<0.98054) {
	agegroups[0]=1;agegroups[1]=1;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.98104) {
	agegroups[0]=1;agegroups[1]=2;agegroups[3]=3;
      } else if (r<0.98223) {
	agegroups[0]=1;agegroups[1]=2;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.98259) {
	agegroups[0]=1;agegroups[1]=2;agegroups[2]=2;agegroups[3]=1;
      } else if (r<0.98668) {
	agegroups[0]=1;agegroups[1]=3;agegroups[3]=2;
      } else if (r<0.98745) {
	agegroups[0]=1;agegroups[1]=3;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.98788) {
	agegroups[0]=2;agegroups[1]=1;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.98987) {
	agegroups[0]=2;agegroups[1]=2;agegroups[3]=2;
      } else if (r<0.99054) {
	agegroups[0]=2;agegroups[1]=2;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.99098) {
	agegroups[0]=2;agegroups[1]=2;agegroups[2]=2;
      } else {
	agegroups[0]=3;agegroups[1]=1;agegroups[3]=2;
      }
    } else { // seven-person household
      if (r<0.99147) {
	agegroups[1]=2;agegroups[2]=2;agegroups[3]=3;
      } else if (r<0.9917) {
	agegroups[1]=2;agegroups[2]=3;agegroups[3]=2;
      } else if (r<0.99185) {
	agegroups[1]=3;agegroups[3]=4;
      } else if (r<0.99205) {
	agegroups[1]=3;agegroups[2]=1;agegroups[3]=3;
      } else if (r<0.99255) {
	agegroups[1]=3;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.99272) {
	agegroups[1]=4;agegroups[3]=2;agegroups[4]=1;
      } else if (r<0.99298) {
	agegroups[1]=4;agegroups[3]=3;
      } else if (r<0.99366) {
	agegroups[1]=4;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.99514) {
	agegroups[1]=5;agegroups[3]=2;
      } else if (r<0.9953) {
	agegroups[1]=5;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.99552) {
	agegroups[0]=1;agegroups[1]=1;agegroups[2]=3;agegroups[3]=2;
      } else if (r<0.99567) {
	agegroups[0]=1;agegroups[1]=2;agegroups[2]=1;agegroups[3]=3;
      } else if (r<0.9961) {
	agegroups[0]=1;agegroups[1]=2;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.9963) {
	agegroups[0]=1;agegroups[1]=3;agegroups[3]=3;
      } else if (r<0.9968) {
	agegroups[0]=1;agegroups[1]=3;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.99694) {
	agegroups[0]=1;agegroups[1]=3;agegroups[2]=2;agegroups[3]=1;
      } else if (r<0.99814) {
	agegroups[0]=1;agegroups[1]=4;agegroups[3]=2;
      } else if (r<0.99836) {
	agegroups[0]=1;agegroups[1]=4;agegroups[2]=1;agegroups[3]=1;
      } else if (r<0.99854) {
	agegroups[0]=2;agegroups[1]=1;agegroups[2]=2;agegroups[3]=2;
      } else if (r<0.99874) {
	agegroups[0]=2;agegroups[1]=2;agegroups[2]=1;agegroups[3]=2;
      } else if (r<0.99889) {
	agegroups[0]=2;agegroups[1]=2;agegroups[2]=2;agegroups[3]=1;
      } else if (r<0.99964) {
	agegroups[0]=2;agegroups[1]=3;agegroups[3]=2;
      } else if (r<0.99984) {
	agegroups[0]=2;agegroups[1]=3;agegroups[2]=1;agegroups[3]=1;
      } else {
	agegroups[0]=3;agegroups[1]=2;agegroups[3]=2;
      }
    }

    int family_size = agegroups[0] + agegroups[1] + agegroups[2] + agegroups[3] + agegroups[4];
    int nNeighborhood = get_rand_uint32 % 4; // random int from 0 to 3
    for (int age=0; age<TAG; age++) {
      while (agegroups[age]>0) {
	int nSchoolgroup = -1;
	if (age==0) { // assign preschool or playgroup?
	  if (get_rand_double<14.0/34.0)  // 14/34 probability of daycare
	    nSchoolgroup = 5+nNeighborhood;
	  else { // otherwise playgroups, 4 children each
	    if (nplay[nNeighborhood] == 4) {  // full, start a new one
	      playgroup[nNeighborhood] = nPlaygroupIndex++;
	      nplay[nNeighborhood] = 0;
	    }
	    nSchoolgroup = playgroup[nNeighborhood];
	    nplay[nNeighborhood]++;
	  }
	} else if (age==1) { // assign school
	  double r2 = get_rand_double;
	  if (r2<0.36)
	    nSchoolgroup = 3 + (nNeighborhood / 2);  // elementary school
	  else if (r2<0.68)
	    nSchoolgroup = 2;  // middle school
	  else if (r2<0.93)
	    nSchoolgroup = 1;  // high school
	  else
	    nSchoolgroup = 0;  // not in school, presumably 18-year-olds or some home-schooled
	}
	create_person(age, family_size, nNumFamilies, nNumFamilies/Community::FAMILIESPERCLUSTER, nNeighborhood, nSchoolgroup, comm);
	comm.nNumResidents++;
	comm.nNumAge[age]++;
	agegroups[age]--;
      }
    }
    nNumFamilies++;
  }
  comm.nLastPerson=nNumPerson;
}

/*
 * Create virtual population
 */
void EpiModel::initPopulation(void) {
  // read census tracts data
  read_tracts();
  // create workgroups
  read_workflow();
  // rescale the contact probabilities
  vector< Community >::iterator cend = commvec.end();
  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    double nightscale = Community::TARGETCOMMUNITYSIZE/(double)comm.nNumResidents; // nighttime population
    double dayscale = Community::TARGETCOMMUNITYSIZE/((double)comm.nNumNonWorkers+comm.nNumWorkers); // daytime population
    for (int i=0; i<5; i++) {
      comm.cpcm[i] = cpcm[i]*nightscale;
      comm.cpnh[i] = cpnh[i]*nightscale;
      comm.daycpcm[i] = cpcm[i]*dayscale;
      comm.daycpnh[i] = cpnh[i]*dayscale;
    }
    for (int i=0; i<9; i++)
      comm.cps[i] = cps[i]*nightscale; // school scaling is based on residential pop
    comm.cps[9] = cps[9];  // note that playgroups are not rescaled
  }
}

/*
 * isAscertained - has this person been ascertained as ill?
 */
bool EpiModel::isAscertained(const Person &p) {
  return isInfectious(p) && getWillBeSymptomatic(p) && getWillBeAscertained(p) && p.iday>=getIncubationDays(p)+nAscertainmentDelay;
}

/*
 * infect
 * set status of "p" to infected
 */
void EpiModel::infect(Person& p) {
  clearSusceptible(p); // no longer susceptible
  setInfected(p);           // infected
  // length of incubation period
  p.iday=-1;  // set to -1 so the person is not infectious until tomorrow
  p.ibits = 0;
  double fSymptomaticProb=fBaseSymptomaticProb;
  if (isVaccinated(p)) {
    if (needsBoost(p))
      fSymptomaticProb*=(1.0-VaccineData[whichVaccine(p)].VEp*defaultvacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
    else 
      fSymptomaticProb*=(1.0-VaccineData[whichVaccine(p)].VEp*VaccineData[whichVaccine(p)].vacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
  }
  if (isAntiviral(p))
    fSymptomaticProb*=(1.0-AVEp);
  double rn = get_rand_double;
  if (rn<fSymptomaticProb) {  // will be symptomatic
    setWillBeSymptomatic(p);
    double rn2 = get_rand_double;
    for (int i=0; i<INCUBATIONMAX; i++) {
      if (rn2 < incubationcdf[i]) {
        setIncubationDays(p,i+1);
	break;
      }
    }
    assert(getIncubationDays(p)>0);
    if (rn<fSymptomaticProb*fSymptomaticAscertainment) { // will be ascertained
      setWillBeAscertained(p);
    }
    rn2 = get_rand_double;
    double *wdcdf = withdrawcdf[isChild(p)?p.age:2];
    if (rn2<wdcdf[0])
      setWithdrawDays(p,getIncubationDays(p)+0);
    else if (rn2<wdcdf[1])
      setWithdrawDays(p,getIncubationDays(p)+1);
    else if (rn2<wdcdf[2])
      setWithdrawDays(p,getIncubationDays(p)+2);
    else
      setWithdrawDays(p,0); // will not withdraw

    if (bTrigger && nTimer>=nTriggerTime &&
	(getWithdrawDays(p)==0 || // doesn't voluntarily withdraw
	 getWithdrawDays(p)-getIncubationDays(p)>1)) { // would withdraw after more than one day
      if ((fLiberalLeaveCompliance>0.0 &&
           isLiberalLeave(tractvec[p.nDayTract-nFirstTract]) && // liberal leave in effect here
	   isWorkingAge(p) && p.nWorkplace>0 && get_rand_double<fLiberalLeaveCompliance) || // will take liberal leave
	  (fVoluntaryIsolationCompliance>0.0 && get_rand_double<fVoluntaryIsolationCompliance)) { // voluntary isolation (not ascertained isolation)
	setWithdrawDays(p,getIncubationDays(p)+1); // stay home the day after symptom onset
      }
    }
    assert(getWillBeSymptomatic(p));
  } else {                // will NOT be symptomatic
    setIncubationDays(p,0);
    setWithdrawDays(p,0);   // note: withdraw days is only checked when symptomatic,
                            // so 0 withdraw days should be ok.
  }
  if (p.nTravelTimer<=0) {
    {
      ++commvec[p.nHomeComm].ninf[p.age];
      ++commvec[p.nHomeComm].nEverInfected[p.age];
    }
  }
}

/*
 * infect
 * "source" infects "p" with contact probability "baseprob"
 * assumes that source is INFECTIOUS
 * returns true if infection occurs, false otherwise
 */
bool EpiModel::infect(Person& p, const Person& source, double baseprob, int sourcetype) {
  assert(source.iday>=0);
  assert(isSusceptible(p));
  if (get_rand_double<baseprob*p.prs*source.pri) {
    p.sourceid = source.id;
    p.sourcetype = sourcetype;
    p.nInfectedTime = nTimer;
    infect(p);
    return true;
  }
  return false;
}

/*
 * isEligible
 * is person p allowed to take vaccine nVacNum
 */
bool EpiModel::isEligible(const Person &p, int nVacNum) {
  if (VaccineData[nVacNum].fAge[p.age]>=1.0 ||
      (VaccineData[nVacNum].bPregnant && isPregnant(p)) ||
      (isInfant(p) && 
       (VaccineData[nVacNum].fInfants>=1.0 ||
	get_rand_double<VaccineData[nVacNum].fInfants)) ||
      get_rand_double<VaccineData[nVacNum].fAge[p.age])
    return false;
  return true;
}

/*
 * vaccinate
 * is person p eligible to get vaccine?
 */
void EpiModel::vaccinate(Person& p) {
  if (!isBoosted(p) && p.nVaccinePriority>0 && !isAscertained(p) && p.nTravelTimer<=0) {
    p.bWantVac=true;
    nNumWantVaccine++;
    assert(p.nTravelTimer<=0);
  }
}

/*
 * resetAscertained
 * reset the number of "ascertained" people in tract t to 0
 */
void EpiModel::resetAscertained(Tract& t, int agegroup=-1) {
  for (unsigned int commid=t.nFirstCommunity; commid<t.nLastCommunity; commid++) {
    for (unsigned int commid=t.nFirstCommunity; commid<t.nLastCommunity; commid++) {
      Community &comm = commvec[commid];
      if (agegroup>=0)
	comm.nRecentlyAscertained[agegroup] = 0;
      else
	memset(comm.nRecentlyAscertained, 0, TAG*sizeof(int));
    }
  }
}

/*
 * vaccinate
 * vaccinate proportion of tract t
 */
void EpiModel::vaccinate(Tract& t) {
  if (!isVaccinated(t)) {
    setVaccinated(t);
    for (unsigned int commid=t.nFirstCommunity; commid<t.nLastCommunity; commid++) {
      Community &comm = commvec[commid];
      for (unsigned int pid=comm.nFirstPerson;
	   pid<comm.nLastPerson;
	   pid++) {
	Person &p = pvec[pid];
	  vaccinate(p);
      }
    }
  }
}

/*
 * TAP
 * TAP person p (antivirals to person and family)
 */
void EpiModel::TAP(Person& p) {
  Community c = commvec[p.nDayComm];
  if (isAVProphylaxis(p))
    clearAVProphylaxis(p); // this person was preparing to get prophylaxis, but now should get treatment
  enum antiviralPolicy policy = getAVPolicy(tractvec[c.nTractID-nFirstTract]);
  if (policy==HHTAP || policy==FULLTAP || (policy==HHTAP100&&nNumTAPDone<100)) {
    nNumTAPDone++;
    for (unsigned int i=c.nFirstPerson; i<c.nLastPerson; i++) {
      Person &target = pvec[i];
      if (!target.bWantAV && 
	  !isAntiviral(target) &&       // not already on antivirals
	  !isAVProphylaxis(target) &&   // not already waiting for AV prophylaxis
	  (p.family == target.family || // family member
	   (policy==FULLTAP && 
	    (p.householdcluster == target.householdcluster || // in household cluster
	     (isChild(p) && isChild(target) && p.nWorkplace==target.nWorkplace && 
	      (p.nWorkplace>=5 ||  // in daycare or playgroup
	       get_rand_double<fContactAscertainment)))))) { // in elementary/middle/high
	if (p.id!=target.id)
	  setAVProphylaxis(target); // prophylaxis, not treatment
	// 100% of family, household cluster, pre-school group, play group
	// fContactAscertainment of other school groups
	target.bWantAV=true;
	nNumWantAV++;
      }
    }
  } else {
    if (!p.bWantAV && !isAntiviral(p) && !isAVProphylaxis(p)) {
      p.bWantAV=true;
      nNumWantAV++;
    }
  }
}

// dayinfectsusceptibles - called by "day" for daytime transmission
void EpiModel::dayinfectsusceptibles(const Person &infected, Community &comm) {
  const double *cpf = (isChild(infected)?cpfc:cpfa); // probability of family transmission
  bool bInfectedIsAtHome = (isWithdrawn(infected) ||
			    isQuarantined(infected) ||
			    isWorkingFromHome(infected) ||
			    infected.nWorkplace==0 ||
			    (isChild(infected) && 
			     infected.nWorkplace<9 &&
			     isSchoolClosed(tractvec[infected.nDayTract-nFirstTract], infected.nWorkplace)));   // infected is at home during the day
  bool bInfectedIsAtSchool = (isChild(infected) && 
			      !isWithdrawn(infected) &&
			      !isQuarantined(infected) &&
			      infected.nWorkplace>0 && 
			      ((infected.age==0 && infected.nWorkplace>=9) ||
			       !isSchoolClosed(tractvec[infected.nDayTract-nFirstTract], infected.nWorkplace))); // infected's school or playgroup is open (playgroups are always open)
  bool bInfectedIsAtWork = (isWorkingAge(infected) &&
			    !isWithdrawn(infected) &&
			    !isQuarantined(infected) &&
			    !isWorkingFromHome(infected) &&
			    infected.nWorkplace>0);  // infected works during the day

  vector< unsigned int >::iterator wend = comm.workers.end();
  list< Person >::iterator vend=comm.visitors.end();
  double casualmultiplier = 1.0; // casual contacts multiplier for children
  if (((infected.age==1 || (infected.age==0 && infected.nWorkplace<9)) &&
       isSchoolClosed(tractvec[infected.nDayTract-nFirstTract],infected.nWorkplace)) ||
      (isWorkingFromHome(infected)))
    casualmultiplier = 2.0;      // casual contacts double for out-of-school children and adults working from home

  // loop over susceptible residents
  for (unsigned int pid2=comm.nFirstPerson;
       pid2<comm.nLastPerson;
       pid2++) {
    Person &p2 = pvec[pid2];
    if (isSusceptible(p2) && p2.nDayComm==comm.id && p2.nTravelTimer<=0) {
      if (infected.family==p2.family && // daytime transmission within household
	  bInfectedIsAtHome &&
	  (isQuarantined(p2) ||
	   isWorkingFromHome(p2) ||
	   p2.nWorkplace==0 ||
	   (isChild(p2) && 
	    p2.nWorkplace<9 &&
	    isSchoolClosed(tractvec[comm.nTractID-nFirstTract], p2.nWorkplace)))) {         // susceptible is at home
	if (infect(p2,infected,cpf[p2.age],FROMFAMILY))
	  continue;
      }

      double casualmultiplier2 = casualmultiplier; // casual contacts multiplier when school is cancelled or workplace is closed
      if (casualmultiplier2==1.0 &&
	  (((p2.age==1 || (p2.age==0 && p2.nWorkplace<9)) && isSchoolClosed(tractvec[p2.nDayTract-nFirstTract],p2.nWorkplace)) ||
	   isWorkingFromHome(p2)))
	casualmultiplier2 = 2.0;      // susceptible's casual contacts double for out-of-school children and adults working from home

      if (!isQuarantined(p2) && !infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.daycpcm[p2.age]*casualmultiplier2,FROMCOMMUNITY)) { // transmission within community
	if (infected.nDayNeighborhood==p2.nDayNeighborhood)  // transmission within neighborhood
	  if (infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.daycpnh[p2.age]*casualmultiplier2,FROMNEIGHBORHOOD))
	    continue;
	if (isChild(infected)) {  // transmitter is child
	  if (bInfectedIsAtSchool && isChild(p2) && infected.nWorkplace==p2.nWorkplace)
	    infect(p2,infected,(infected.nWorkplace>=9?cps[9]:comm.cps[infected.nWorkplace]),FROMSCHOOL); // transmission within school/daycare/playgroup
	} else {           // transmitter is adult
	  if (bInfectedIsAtWork && isWorkingAge(p2) && !isWorkingFromHome(p2) && infected.nWorkplace==p2.nWorkplace)  // transmission within work group
	    infect(p2,infected,cpw,FROMWORK);
	}
      }
    }
  }

  // loop over susceptible workers visiting from other communities
  for (vector< unsigned int >::iterator it = comm.workers.begin(); 
       it != wend;
       it++) {
    Person &p2 = pvec[*it];
    if (isSusceptible(p2) && !isQuarantined(p2) && !isWorkingFromHome(p2) && p2.nTravelTimer<=0) {
      if (!infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.daycpcm[p2.age]*casualmultiplier,FROMCOMMUNITY)) { // transmission within community
	if (infected.nDayNeighborhood==p2.nDayNeighborhood)  // transmission within neighborhood
	  if (infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.daycpnh[p2.age]*casualmultiplier,FROMNEIGHBORHOOD))
	    continue;
	if (isWorkingAge(infected) && infected.nWorkplace==p2.nWorkplace) {
	  // transmit to coworkers from other tracts
	  assert(infected.nWorkplace>0);
	  infect(p2,infected,cpw,FROMWORK);
	}
      }
    }
  }

  // loop over susceptible visitors (short term travelers)
  if (bTravel) {
    for (list< Person >::iterator it = comm.visitors.begin();
	 it != vend;
	 it++) {
      Person &p2 = *it;
      assert(p2.nDayComm==comm.id);
      assert(p2.nTravelTimer>0);
      if (isSusceptible(p2)) { // && !isQuarantined(p2)) {
	if (!infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.daycpcm[p2.age]*casualmultiplier,FROMCOMMUNITY)) { // transmission within community
	  if (infected.nDayNeighborhood==p2.nDayNeighborhood)  // transmission work neighborhood
	    if (infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.daycpnh[p2.age]*casualmultiplier,FROMNEIGHBORHOOD))
	      continue;
	  if (isWorkingAge(infected) && isWorkingAge(p2) && infected.nWorkplace==p2.nWorkplace && p2.nWorkplace>0)
	    infect(p2,infected,cpw,FROMWORK);
	}
      }
    }
  }
}

/*
 * day
 * daytime transmission
 */
void EpiModel::day(void) {
  // calculate daytime susceptibility and infectiousness for each person
  vector< Community >::iterator cend = commvec.end();
  vector< Person >::iterator pend=pvec.end();
  for (vector< Person >::iterator it = pvec.begin(); 
       it != pend;
       it++) {
    Person &p = *it;
    p.prs = 1.0-p.fBaselineVES; // susceptibility multiplier
    if (isInfectious(p))
      p.pri = vload[p.nWhichVload][(int)(p.iday)] *
	seasonality[nTimer/2]; // infectiousness multiplier
    else
      p.pri = 0.0;
    if (isVaccinated(p)) {  // effects of the vaccine on transmission
      if (needsBoost(p)) {
	p.prs *= (1.0-VaccineData[whichVaccine(p)].VEs*defaultvacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
	p.pri *= (1.0-VaccineData[whichVaccine(p)].VEi*defaultvacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
      } else {
	p.prs *= (1.0-VaccineData[whichVaccine(p)].VEs*VaccineData[whichVaccine(p)].vacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
	p.pri *= (1.0-VaccineData[whichVaccine(p)].VEi*VaccineData[whichVaccine(p)].vacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
      }
    }
    if (isAntiviral(p)) {  // effects of the antiviral on transmission
      p.prs *= (1.0-AVEs);
      p.pri *= (1.0-AVEi);
    }
    if (isSymptomatic(p))
      p.pri *= fRelativeSymptomaticInfectiousness;  // symptomatic people are 2 times as infectious
    assert(p.pri<=1.0);
    assert(p.prs<=1.0);
  }

  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    list< Person >::iterator vend=comm.visitors.end();
    for (list< Person >::iterator pit = comm.visitors.begin(); 
	 pit != vend;
	 pit++) {
      Person &p = *pit;
      p.prs = 1.0-p.fBaselineVES; // susceptibility multiplier
      if (isInfectious(p))
	p.pri = vload[p.nWhichVload][(int)(p.iday)]; // infectiousness multiplier
      else
	p.pri = 0.0;
      if (isVaccinated(p)) {
	if (needsBoost(p)) {
	  p.prs *= (1.0-VaccineData[whichVaccine(p)].VEs*defaultvacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
	  p.pri *= (1.0-VaccineData[whichVaccine(p)].VEi*defaultvacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
	} else {
	  p.prs *= (1.0-VaccineData[whichVaccine(p)].VEs*VaccineData[whichVaccine(p)].vacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
	  p.pri *= (1.0-VaccineData[whichVaccine(p)].VEi*VaccineData[whichVaccine(p)].vacceff[p.vday]*fVaccineEfficacyByAge[p.age]);
	}
      }
      if (isAntiviral(p)) {  // effects of the antiviral on transmission
	p.prs *= (1.0-AVEs);
	p.pri *= (1.0-AVEi);
      }
      if (isSymptomatic(p))
	p.pri *= fRelativeSymptomaticInfectiousness;
      assert(p.pri<=1.0);
      assert(p.prs<=1.0);
    }
  }

  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    if (comm.ninf[0]>0 || comm.ninf[1]>0 || comm.ninf[2]>0 || comm.ninf[3]>0 || comm.ninf[4]>0) { // is any resident currently infected?
      // loop over infectious residents
      for (unsigned int pid=comm.nFirstPerson;
	   pid<comm.nLastPerson;
	   pid++) {
	Person &p = pvec[pid];
	if (isInfectious(p) && !isWithdrawn(p) && !isQuarantined(p) && p.nDayComm==comm.id && p.nTravelTimer<=0) {
	  dayinfectsusceptibles(p, comm);
	}
      }
    }

    // loop over infectious workers visiting from other communities on this processor
    vector< unsigned int >::iterator wend = comm.workers.end();
    for (vector< unsigned int >::iterator it = comm.workers.begin();
	 it != wend;
	 it++) {
      Person &p = pvec[*it];
      if (isInfectious(p) && !isWithdrawn(p) && !isQuarantined(p) && !isWorkingFromHome(p) && p.nTravelTimer<=0)
	dayinfectsusceptibles(p, comm);
    }

    // loop over infectious visitors (short-term travelers)
    if (bTravel) {
      list< Person >::iterator vend = comm.visitors.end();
      for (list< Person >::iterator it = comm.visitors.begin(); 
	   it != vend;
	   it++) {
	Person &p = *it;
	assert(p.nDayComm==comm.id);
	if (isInfectious(p) && !isWithdrawn(p)) // && !isQuarantined(p))
	  dayinfectsusceptibles(p, comm);
      }
    }
  }
}

// nightinfectsusceptibles - called by "night" for nighttimetime transmission
void EpiModel::nightinfectsusceptibles(const Person &infected, Community &comm) {
  const double *cpf = (isChild(infected)?cpfc:cpfa); // probability of family transmission
  const double *cphc = (isChild(infected)?cphcc:cphca); // probability of household cluster transmission
  
  for (unsigned int pid2=comm.nFirstPerson;
       pid2<comm.nLastPerson;
       pid2++) {
    Person &p2 = pvec[pid2];
    if (isSusceptible(p2) && p2.nTravelTimer<=0) { // && p.id!=p2.id) {
      if (infected.family==p2.family)   // transmission within household
	if (infect(p2,infected,cpf[p2.age],FROMFAMILYNIGHT))
	  continue;
      if (!isWithdrawn(infected) && !isQuarantined(infected) && !isQuarantined(p2)) {
	assert(infected.nHomeComm==p2.nHomeComm);
	if (!infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.cpcm[p2.age],FROMCOMMUNITYNIGHT)) { // transmission within home community
	  if (infected.nHomeNeighborhood==p2.nHomeNeighborhood) // transmission within neighborhood
	    if (infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.cpnh[p2.age],FROMNEIGHBORHOODNIGHT))
	      continue;
	  if (infected.householdcluster==p2.householdcluster) // transmission within household cluster
	    infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*cphc[p2.age],FROMHHCLUSTERNIGHT);
	}
      }
    }
  }

  // check for susceptible travelers
  if (bTravel) {
    list< Person >::iterator vend=comm.visitors.end();
    for (list< Person >::iterator it = comm.visitors.begin(); 
	 it != vend;
	 it++) {
      Person &p2 = *it;
      if (isSusceptible(p2)) {
	if (infected.family==p2.family)  // transmission within household
	  if (infect(p2,infected,cpf[p2.age],FROMFAMILY))
	    continue;
	if (!isWithdrawn(infected) && !isQuarantined(infected) && !isQuarantined(p2)) {
	  assert(infected.nHomeComm==p2.nHomeComm);
	  if (!infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.cpcm[p2.age],FROMCOMMUNITYNIGHT)) { // transmission within home community
	    if (infected.nHomeNeighborhood==p2.nHomeNeighborhood) // transmission within neighborhood
	      if (infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*comm.cpnh[p2.age],FROMNEIGHBORHOODNIGHT))
		continue;
	    if (infected.householdcluster==p2.householdcluster) // transmission within household cluster
	      infect(p2,infected,(isCommunityContactReduction(tractvec[comm.nTractID-nFirstTract])?(1.0-fCommunityContactReduction):1.0)*cphc[p2.age],FROMHHCLUSTERNIGHT);
	  }
	}
      }
    }
  }
}

/*
 * night
 * nighttime transmission and updating clocks
 */
void EpiModel::night(void) {
  vector< Community >::iterator cend = commvec.end();
  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    if (comm.ninf[0]>0 || comm.ninf[1]>0 || comm.ninf[2]>0 || comm.ninf[3]>0 || comm.ninf[4]>0) { // is any resident currently infected?
      // transmit infection
      for (unsigned int pid=comm.nFirstPerson;
	   pid<comm.nLastPerson;
	   pid++) {
	Person &p = pvec[pid];
	if (isInfectious(p) && p.nTravelTimer<=0)
	  nightinfectsusceptibles(p, comm);
      }
    }

    // check for infectious travelers
    if (bTravel) {
      list< Person >::iterator vend=comm.visitors.end();
      for (list< Person >::iterator it = comm.visitors.begin(); 
	   it != vend;
	   it++) {
	Person &p = *it;
	assert(p.nDayComm==comm.id);
	if (isInfectious(p) && !isWithdrawn(p))
	  nightinfectsusceptibles(p, comm);
      }
    }

    // update infection and intervention status (timers)
    for (unsigned int pid=comm.nFirstPerson;
	 pid<comm.nLastPerson;
	 pid++) {
      Person &p = pvec[pid];
      if (isInfected(p)) {
	if (isInfectious(p)) {
	  if (getWillBeSymptomatic(p) && p.iday==getIncubationDays(p)) {
	    setSymptomatic(p);
	    ++commvec[p.nHomeComm].nEverSymptomatic[p.age];
	    ++commvec[p.nHomeComm].nsym[p.age];
	  }
	  if (isSymptomatic(p)) {
	    if (p.iday==getIncubationDays(p)+nAscertainmentDelay && getWillBeAscertained(p)) { // person becomes ascertained
	      ++commvec[p.nHomeComm].nEverAscertained[p.age];
	      ++commvec[p.nHomeComm].nRecentlyAscertained[p.age];
	      // call TAP when cases are ascertained in this tract
	      if (getAVPolicy(tractvec[comm.nTractID-nFirstTract])!=NOAV)
		TAP(p);
	      if (isAscertainedIsolation(tractvec[comm.nTractID-nFirstTract]) && get_rand_double<fAscertainedIsolationCompliance) { // ascertained isolation
		setWithdrawDays(p,getIncubationDays(p)+nAscertainmentDelay);
	      }
	      
	      // quarantine the family
	      if (isQuarantine(tractvec[comm.nTractID-nFirstTract])) {
		for (unsigned int pid2=comm.nFirstPerson;
		     pid2<comm.nLastPerson;
		     pid2++) {
		  Person &p2 = pvec[pid2];
		  if (p.family==p2.family && p.id!=p2.id && get_rand_double<fQuarantineCompliance) { // quarantine family member
		    setQuarantined(p2);  // household quarantine
		    p2.nQuarantineTimer = nQuarantineLength+1;
		  }
		}
	      }
	    }
	    if (getWithdrawDays(p)>0 && p.iday==(int)getWithdrawDays(p))
	      setWithdrawn(p);              // withdraw to home
	  }
	}
	p.iday++;
	if (p.iday>=VLOADNDAY) {
	  if (isSymptomatic(p))
	    comm.nsym[p.age]--;
	  clearSymptomatic(p); // recovered
	  clearSusceptible(p); // recovered
	  clearInfected(p);    // recovered
	  clearWithdrawn(p);   // recovered
	  comm.ninf[p.age]--;
	  p.iday=0;
	}
      }
      if (isVaccinated(p) && p.vday<VACCEFFLENGTH && 
	  (isBoosted(p) ||
	   p.vday<VaccineData[whichVaccine(p)].nBoostDay ||
	   (!needsBoost(p) &&
	    VaccineData[whichVaccine(p)].nNumDoses==1) ||
	   (needsBoost(p) &&
	    p.vday<defaultboostday)))
	p.vday++;
      if (isQuarantined(p) && --p.nQuarantineTimer<=0) {
	clearQuarantined(p);  // quarantine over for this person
      }
      if (isAntiviral(p)) {
	--p.nAVTimer;
	if (!isAVProphylaxis(p))
	  --p.nAVTimer; // treatment requires 2nd tablet each day
	if (p.nAVTimer<=0 || 
	    (p.nAVTimer==nAntiviralCourseSize-2 && get_rand_double<fStopAntiviralTwoPills)) {
	  clearAntiviral(p);    // antiviral over for this person
	  clearAVProphylaxis(p);
	}
      }
    }

    // update visitors
    if (bTravel) {
      list< Person >::iterator vend=comm.visitors.end();
      for (list< Person >::iterator vit = comm.visitors.begin(); 
	   vit != vend;
	   vit++) {
	Person &p = *vit;
	if (isInfected(p)) {
	  if (isInfectious(p)) {
	    if (p.iday==getIncubationDays(p) && getWillBeSymptomatic(p))
	      setSymptomatic(p);
	    if (getWithdrawDays(p)>0 && p.iday==(int)getWithdrawDays(p))
	      setWithdrawn(p); 	// withdraw to home
	  }
	  p.iday++;
	  if (p.iday>=VLOADNDAY) {
	    clearSymptomatic(p); // recovered
	    clearSusceptible(p); // recovered
	    clearInfected(p);    // recovered
	    clearWithdrawn(p);   // recovered
	    p.iday=0;
	  }
	}
	if (isVaccinated(p) && p.vday<VACCEFFLENGTH && 
	    (VaccineData[whichVaccine(p)].nNumDoses==1 || 
	     p.vday<VaccineData[whichVaccine(p)].nBoostDay || 
	     isBoosted(p)))
	  p.vday++;
      }
    }
  }
}

/*
 * travel_end
 * return from short-term travel to a random census tract
 */
void EpiModel::travel_end(void) {
  // update travel timers and send some people home
  vector< Community >::iterator cend = commvec.end();
  for (vector< Community >::iterator it = commvec.begin(); 
       it != cend;
       it++) {
    Community &comm = *it;
    list< Person >::iterator vend=comm.visitors.end();
    for (list< Person >::iterator it = comm.visitors.begin(); 
	 it != vend;
	 ) {
      Person &p = *it;
      if (--p.nTravelTimer<=0) {
	{
	  // update infection status in the returning person
	  assert(p.id<pvec.size());
	  Person &original = pvec[p.id];
	  if (isInfected(p)) { // got sick on vacation
	    if (!isInfected(original)) {
	      ++commvec[original.nHomeComm].ninf[original.age];
	      ++commvec[original.nHomeComm].nEverInfected[original.age];
	      original.nInfectedTime = p.nInfectedTime;
	    }
	    if (isSymptomatic(p) && !isSymptomatic(original)) {
	      ++commvec[original.nHomeComm].nEverSymptomatic[original.age];
	      ++commvec[original.nHomeComm].nsym[original.age];
	    }
	  } else if (isInfected(original)) {
	    --commvec[original.nHomeComm].ninf[original.age]; // got better on vacation
	    if (!isSymptomatic(p) && isSymptomatic(original))
	      --commvec[original.nHomeComm].nsym[original.age];
	  }
	  original.status = p.status;
	  original.ibits = p.ibits;
	  original.iday = p.iday;
	  original.sourcetype = p.sourcetype;
	  original.sourceid = p.sourceid;
	  original.nTravelTimer = 0; // welcome home
	}
	it = comm.visitors.erase(it);
	vend = comm.visitors.end();
      } else
	it++;
    }
  }
}

/*
 * travel_start
 * send people on short-term travel to a random census tract
 */
void EpiModel::travel_start(void) {
  // update travel timers and send some people home
  // send people on trips
  // Quarantined people don't travel.
  // Recovered/immune people probably shouldn't travel.
  vector< Person >::iterator pend = pvec.end();
  for (vector< Person >::iterator it = pvec.begin(); 
       it != pend;
       it++) {
    Person &p = *it;
    // This is inefficient!  try using a binomial to compute the number
    // of travelers.
    if (p.nTravelTimer<=0 && !isQuarantined(p) && get_rand_double < travel_pr[p.age]) {
      // We're going to DisneyWorld!
      double r = get_rand_double;
      int i;
      for (i=0; r>travel_length_cdf[i]; i++)
	;
      p.nTravelTimer = i+1;
      unsigned int destinationtract = get_rand_uint32 % nNumTractsTotal;
      {
	// assign traveler to a community in tract
	unsigned int destinationcomm = tractvec[destinationtract-nFirstTract].nFirstCommunity;
	int diff = tractvec[destinationtract-nFirstTract].nLastCommunity-tractvec[destinationtract-nFirstTract].nFirstCommunity;
	if (diff>1)
	  destinationcomm += get_rand_uint32 % diff;
	assert(destinationcomm<commvec.size());
	commvec[destinationcomm].visitors.push_back(p);
	Person &traveler = commvec[destinationcomm].visitors.back();
	// assign new neighborhood, workplace, etc. to the traveler
	traveler.nDayTract = destinationtract;
	traveler.nHomeComm = destinationcomm;
	traveler.nDayComm = destinationcomm;
	traveler.nWorkNeighborhood = get_rand_uint32 % 4; // work neighborhood
	traveler.nDayNeighborhood = traveler.nWorkNeighborhood;
	diff = commvec[destinationcomm].nLastPerson-commvec[destinationcomm].nFirstPerson;
	if (diff>0) {
	  // assign a host family to the traveler
	  Person &host = pvec[commvec[destinationcomm].nFirstPerson + (get_rand_uint32 % diff)];
	  traveler.family = host.family;
	  traveler.householdcluster = host.householdcluster;
	  traveler.nHomeNeighborhood = host.nHomeNeighborhood;
	  if (isWorkingAge(traveler) && commvec[destinationcomm].nNumWorkGroups>0)
	    traveler.nWorkplace = 1 + (get_rand_uint32 % commvec[destinationcomm].nNumWorkGroups);
	  else
	    traveler.nWorkplace = 0;
	}
      }
    }
  }
}

/*
 * response
 * Determine if the epidemic has been detected and
 * make reactive responses if appropriate
 */
void EpiModel::response(void) {
  // 0. Open schools as appropriate
  for (vector< Tract >::iterator it = tractvec.begin();
       it != tractvec.end();
       it++) {
    Tract &t = *it;
    if (isSchoolClosed(t,1) && nSchoolOpeningDays[t.fips_state-1]-1==nTimer/2) {
      for (int i=0; i<9; i++) {
	setSchoolOpen(t,i);// schools open today
      }
      t.nSchoolClosureTimer=0;
      resetAscertained(t,1); // reset school-aged children
    }
  }

  // 1. Count cumulative number of ascertained cases
  if (!bTrigger) {
    int nNumRecentlyAscertained=0; // number of people recently ascertained
    vector< Community >::iterator cend = commvec.end();
    for (vector< Community >::iterator it = commvec.begin();
	 it != cend;
	 it++) {
      Community &c = *it;
      for (int j=0; j<TAG; j++)
	nNumRecentlyAscertained += c.nRecentlyAscertained[j];
    }
    if (nNumRecentlyAscertained/(double)nNumPerson>fResponseThreshold) { // trigger reached
      bTrigger=true;
      nTriggerTime=nTimer+nTriggerDelay*2;
      cout << "Response starts on day " << (nTriggerTime/2) << ": ascertained: " << nNumRecentlyAscertained << "/" << nNumPerson << "=" << nNumRecentlyAscertained/(double)nNumPerson << " on day " << (nTimer/2) << endl;
      for (vector< Tract >::iterator it = tractvec.begin();
	   it != tractvec.end();
	   it++) {
	Tract &t = *it;
	resetAscertained(t);
      }
    }
  }

  // 2. Epidemic is detected, start reactive strategies
  if (bTrigger && nTimer==nTriggerTime) {
    for (vector< Tract >::iterator it = tractvec.begin();
	 it != tractvec.end();
	 it++) {
      Tract &t = *it;
      if (eAVPolicy!=NOAV && getAVPolicy(t)!=eAVPolicy)
	setAVPolicy(t, eAVPolicy);         // start using antivirals in tract
      if (fQuarantineCompliance>0.0 && !isQuarantine(t))
	setQuarantine(t);  // activate household quarantine in tract
      if (fAscertainedIsolationCompliance>0.0 && !isAscertainedIsolation(t))
	setAscertainedIsolation(t);  // activate ascertained isolation in tract
      if (nSchoolClosureDays>0 && nSchoolClosurePolicy==1) {
	cout << "Closing schools in tract " << t.id << " on day " << (nTimer/2) << " for " << nSchoolClosureDays << " days" << endl;
	t.nSchoolClosureTimer=nSchoolClosureDays;
	if (!isSchoolClosed(t,1)) {
	  for (int i=0; i<9; i++)
	    setSchoolClosed(t,i);// activate school closures
	}
      }
      if (fWorkFromHomeCompliance>0.0 && nWorkFromHomeDuration>0 && !isWorkFromHome(t)) {
	setWorkFromHome(t);  // activate work from home in tract
	t.nWorkFromHomeTimer = nWorkFromHomeDuration;
	cout << "Working from home in tract " << t.id << " on day " << (nTimer/2) << " for " << nWorkFromHomeDuration << " days" << endl;
      }
      if (fLiberalLeaveCompliance>0.0 && nLiberalLeaveDuration>0 && !isLiberalLeave(t)) {
	setLiberalLeave(t);  // activate liberal leave in tract
	t.nLiberalLeaveTimer = nLiberalLeaveDuration;
	cout << "Liberal leave in tract " << t.id << " on day " << (nTimer/2) << " for " << nLiberalLeaveDuration << " days" << endl;
      }
      if (fCommunityContactReduction>0.0 && nCommunityContactReductionDuration>0 && !isCommunityContactReduction(t)) {
	setCommunityContactReduction(t);  // activate community contact reduction in tract
	t.nCommunityContactReductionTimer = nCommunityContactReductionDuration;
	cout << "Community contact reduction of " << 100*fCommunityContactReduction << "% in tract " << t.id << " on day " << (nTimer/2) << " for " << nCommunityContactReductionDuration << " days" << endl;
      }
      if (eVaccinationStrategy==MASSVAC && !isVaccinated(t))
	vaccinate(t);
    }

    if (fWorkFromHomeCompliance>0.0) { // choose which people will start working from home
      vector< Person >::iterator pend=pvec.end();
      for (vector< Person >::iterator it = pvec.begin(); 
	   it != pend;
	   it++) {
	Person &p = *it;
	if ((isWorkingAge(p)) && (p.nWorkplace>=0) &&
	    (get_rand_double < fWorkFromHomeCompliance)) {
	  setWorkingFromHome(p);
	  p.nDayNeighborhood = p.nHomeNeighborhood; // now spends daytime at home
	}
      }
    }
  } else {
    assert(nTriggerTime % 2==1); // make sure nTriggerTime is an odd number or else it might not work
  }

  // 3. change vaccine priorities if needed
  if (nPriorityChangeTime>0 && nPriorityChangeTime==nTimer/2) {
    vector< Person >::iterator pend=pvec.end();
    for (vector< Person >::iterator it = pvec.begin(); 
	 it != pend;
	 it++) {
      Person &p = *it;
      if (p.nVaccinePriority>0) { // is this person already able to get vaccine?
	unsigned char priority = 100;
	if (nVaccinePriorities2[PRIORITY_PREGNANT]>0 && isPregnant(p))
	  priority = min(nVaccinePriorities2[PRIORITY_PREGNANT], priority);
	if (nVaccinePriorities2[PRIORITY_ESSENTIAL]>0 && isEssential(p))
	  priority = min(nVaccinePriorities2[PRIORITY_ESSENTIAL], priority);
	if (nVaccinePriorities2[PRIORITY_HR0+p.age]>0 && isHighRisk(p))
	  priority = min(nVaccinePriorities2[PRIORITY_HR0+p.age], priority);
	if (nVaccinePriorities2[PRIORITY_0+p.age]>0)
	  priority = min(nVaccinePriorities2[PRIORITY_0+p.age], priority);
	if (nVaccinePriorities2[PRIORITY_INFANTFAMILY]>0 && !isInfant(p)) {
	  unsigned int upper = (p.id+7>pvec.size()?pvec.size():p.id+7);
	  for (unsigned int i=(p.id>7?p.id-7:0); i<upper; i++) {
	    if (pvec[i].family==p.family &&
		isInfant(pvec[i])) {
	      priority = min(nVaccinePriorities2[PRIORITY_INFANTFAMILY], priority);
	      i=upper;
	    }
	  }
	}
	if (priority>nHighestPriority)
	  p.nVaccinePriority = 0;
	else
	  p.nVaccinePriority = priority;
      }
    }
  }

  // 4. Ring vaccination and consider closing or re-opening schools
  if (bTrigger && nTimer>=nTriggerTime) {
    if (eVaccinationStrategy==RINGVACTRACT) {
      for (vector< Tract >::iterator it = tractvec.begin();
	   it != tractvec.end();
	   it++) {
	Tract &t = *it;
	if (!isVaccinated(t)) {
	  int nNumRecentlyAscertained=0; // number of people ascertained in this tract
	  for (unsigned int i=t.nFirstCommunity; i<t.nLastCommunity; i++)
	    for (int j=0; j<TAG; j++)
	      nNumRecentlyAscertained += commvec[i].nRecentlyAscertained[j];
	  if (nNumRecentlyAscertained>0) { // case found in tract
	    vaccinate(t);
	  }
	}
      }
    }
    if (nSchoolClosureDays>0) { // is school closure an option?
      if (nSchoolClosurePolicy==2) { // close schools by ascertainment
	for (vector< Tract >::iterator it = tractvec.begin();
	     it != tractvec.end();
	     it++) {
	  Tract &t = *it;
	  int nNumRecentlyAscertained=0; // number of children ascertained in this tract
	  for (unsigned int i=t.nFirstCommunity; i<t.nLastCommunity; i++)
	    nNumRecentlyAscertained += commvec[i].nRecentlyAscertained[1]; // count school-aged children only
	  if (nNumRecentlyAscertained>0 && 
	      (!isSchoolClosed(t,1) || !isSchoolClosed(t,2) ||
	       !isSchoolClosed(t,3) || !isSchoolClosed(t,4) ||
	       !isSchoolClosed(t,5) || !isSchoolClosed(t,6) ||
	       !isSchoolClosed(t,7) || !isSchoolClosed(t,8))) { // close schools in this tract if any are still open
	    t.nSchoolClosureTimer=nSchoolClosureDays;
	    cout << "Closing schools in tract " << t.id << " on day " << (nTimer/2) << endl;
	    for (int i=0; i<9; i++) // close all schools in the tract
	      setSchoolClosed(t,i);
	  }
	}
      }

      for (vector< Tract >::iterator it = tractvec.begin();
	   it != tractvec.end();
	   it++) {
	Tract &t = *it;
	if (isSchoolClosed(t,1) &&                       // school is closed
	    --t.nSchoolClosureTimer<=0 &&                // school closure is not in effect anymore
	    nSchoolOpeningDays[t.fips_state-1]-1<=nTimer/2) { // school should be open
	  cout << "School open in tract " << t.id << " on day " << nTimer/2 << endl;
	  for (int i=0; i<9; i++)
	    setSchoolOpen(t,i);// school is open again
	  resetAscertained(t,1); // reset school-aged children
	}
      }
    }
    // are non-school NPIs done?
    if ((fLiberalLeaveCompliance>0.0 && nLiberalLeaveDuration>0 && nLiberalLeaveDuration<5000) ||
	(fWorkFromHomeCompliance>0.0 && nWorkFromHomeDuration>0 && nWorkFromHomeDuration<5000) ||
        (fCommunityContactReduction>0.0 && nCommunityContactReductionDuration>0 && nCommunityContactReductionDuration<5000)) { // is liberal leave or work from home an option and of finite duration?
      bool bBackToWork = false;
      for (vector< Tract >::iterator it = tractvec.begin();
	   it != tractvec.end();
	   it++) {
	Tract &t = *it;
	if (isWorkFromHome(t) &&
	    (--t.nWorkFromHomeTimer)<=0) {
	  //	  cout << "no more working at home in tract " << t.id << " on day " << nTimer/2 << endl;
	  clearWorkFromHome(t);// tract goes back to work
	  bBackToWork = true;
	}
	if (isLiberalLeave(t) &&
	    (--t.nLiberalLeaveTimer)<=0) {
	  //	  cout << "need to work when sick in tract " << t.id << " on day " << nTimer/2 << endl;
	  clearLiberalLeave(t);// sick people in this tract don't go home as much
	}
	if (isCommunityContactReduction(t) &&
	    (--t.nCommunityContactReductionTimer)<=0) {
	  clearCommunityContactReduction(t);// resume community contacts in this tract
	}
      }

      // really inefficient loop to send people back to work
      if (bBackToWork) { // send people back to work
	vector< Person >::iterator pend=pvec.end();
	for (vector< Person >::iterator it = pvec.begin(); 
	     it != pend;
	     it++) {
	  Person &p = *it;
	  if (isWorkingAge(p) && (p.nWorkplace>=0) &&
	      isWorkingFromHome(p) &&
	       !isWorkFromHome(tractvec[p.nDayTract-nFirstTract])) {
	    clearWorkingFromHome(p);
	    p.nDayNeighborhood = p.nWorkNeighborhood; // now spends daytime at work
	  }
	}
      }
    }
  }

  // 5. Count available vaccines
  for (int i=0; i<nNumVaccineTypes; i++)
    if (vaccineproductionschedule[i][nTimer/2]>0)
      nVaccineSupply[i] += vaccineproductionschedule[i][nTimer/2];  // new vaccine arrives

  unsigned int nVacLim[NUMVACCINES];
  unsigned int totalvaccinesupply=0;
  for (int i=0; i<nNumVaccineTypes; i++) {
    nVacLim[i] = nVaccineSupply[i]-nNumVaccineDosesUsed[i]; // how many vaccine doses are available today?
    if (nVaccineSupply[i]<=nNumVaccineDosesUsed[i]) // because we are using unsigned ints
      nVacLim[i]=0;
    totalvaccinesupply += nVacLim[i];
  }
  if (totalvaccinesupply>nVaccineDailyLimit) { // limited overall vaccine daily delivery capacity
    double factor = nVaccineDailyLimit/(double)totalvaccinesupply;
    for (int i=0; i<nNumVaccineTypes; i++)
      nVacLim[i] *= factor;
  }

  // 6. Distribute vaccines
  vector< Person >::iterator pend=pvec.end();
  for (int vacnum=0; vacnum<nNumVaccineTypes; vacnum++) {
    if (nVacLim[vacnum]>0) {
      double vFrac=0.0;
      // demand count for THIS vaccine by priority tier
      unsigned int wantthisvac[PRIORITY_LAST];
      memset(wantthisvac, 0, PRIORITY_LAST*sizeof(unsigned int));
      for (vector< Person >::iterator it = pvec.begin(); 
	   it != pend;
	   it++) {
	Person &p = *it;
	if (p.bWantVac &&
	    p.bVaccineEligible[vacnum] &&
	    (!isVaccinated(p) ||
	     (whichVaccine(p)==vacnum &&
	      ((VaccineData[vacnum].nNumDoses>1 &&
		p.vday>=VaccineData[vacnum].nBoostDay) ||
	       (needsBoost(p) && p.vday>=defaultboostday)))))
	  wantthisvac[p.nVaccinePriority]++;
      }
      for (int priority=1; priority<=nHighestPriority; priority++) {
	if (wantthisvac[priority]>0) {
	  vFrac = nVacLim[vacnum]/(double)wantthisvac[priority];
	  for (vector< Person >::iterator it = pvec.begin(); 
	       it != pvec.end();
	       it++) {
	    Person &p = *it;
	    if (p.bWantVac &&
		getVaccinePriority(p)==priority &&
		p.bVaccineEligible[vacnum] &&
		(!isVaccinated(p) ||
		 (whichVaccine(p)==vacnum &&
		  (p.vday>=VaccineData[vacnum].nBoostDay ||
		   (needsBoost(p) && p.vday>=defaultboostday))))) {
	      if (vFrac>=1.0 || get_rand_double < vFrac) {
		if (needsBoost(p) && VaccineData[vacnum].nNumDoses>1)
		  clearNeedsBoost(p);
		nNumVaccineDosesUsed[vacnum]++;
		if (nVacLim[vacnum]>0)
		  nVacLim[vacnum]--;
		setWhichVaccine(p, vacnum);
		if (!isVaccinated(p)) {
		  setVaccinated(p);
		  p.vday = 0;
		  if (VaccineData[vacnum].nNumDoses<=1 && !needsBoost(p)) {
		    p.bWantVac = false;
		    setBoosted(p);
		    nNumWantVaccine--;
		  } else if (ePrevaccinationStrategy==PRIMEBOOSTRANDOM) {
		    nNumWantVaccine--;
		    p.bWantVac = false; // only gets one shot if PRIMEBOOSTRANDOM
		  }
		} else {
		  p.bWantVac = false;
		  setBoosted(p);
		  nNumWantVaccine--;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // 7. Count available antiviral courses
  unsigned int uAVFraction;
  unsigned int nAVLim = min(nAVTotalLimit-nNumAntiviralsUsed, nAVDailyLimit); // how many antiviral courses are available today?
  if (nAVTotalLimit<=nNumAntiviralsUsed) // because we are using unsigned ints
    nAVLim=0;
  if (nNumWantAV<=nAVLim)
    uAVFraction=UINT_MAX;
  else
    uAVFraction=(unsigned int)(nAVLim/(double)nNumWantAV*UINT_MAX);

  // 7. Distribute antivirals
  for (vector< Person >::iterator it = pvec.begin(); 
       it != pend;
       it++) {
    Person &p = *it;
    if (p.bWantAV && (uAVFraction==UINT_MAX || get_rand_uint32<uAVFraction)) {
      p.bWantAV = false;
      setAntiviral(p);
      p.nAVTimer = nAntiviralCourseSize;  // supply course of AV to this person
      nNumAntiviralsUsed++;
      nNumWantAV--;
    }
  }

  // 8. Record number of vaccine and antivirals distributed
  for (int i=0; i<nNumVaccineTypes; i++)
    nNumVaccineDosesUsedHistory[i].push_back(nNumVaccineDosesUsed[i]);
  nNumAntiviralsUsedHistory.push_back(nNumAntiviralsUsed);
}

// log - outputs data from the current timestep to a log file
void EpiModel::log(void) {
  if (nLogFileInterval<=0)
    return;
  ostream &out = *logfile;
  // output number of infecteds by tract
  for (vector< Tract >::iterator it = tractvec.begin(); 
       it != tractvec.end();
       it++) {
    Tract& t = *it;
    out << int(nTimer/2) << "," << t.id;
    int ninf[TAG],      // current infected prevalence
      ncinf[TAG],	// cumulative infected attack rate
      nsym[TAG],      // current symptomatic prevalence
      ncsym[TAG];	// cumulative symptomatic attack rate
    memset(ninf, 0, sizeof(int)*TAG);
    memset(ncinf, 0, sizeof(int)*TAG);
    memset(nsym, 0, sizeof(int)*TAG);
    memset(ncsym, 0, sizeof(int)*TAG);
    for (unsigned int i=t.nFirstCommunity; i<t.nLastCommunity; i++) {
      for (int j=0; j<TAG; j++) {
	ninf[j] += commvec[i].ninf[j];
	ncinf[j] += commvec[i].nEverInfected[j];
	nsym[j] += commvec[i].nsym[j];
	ncsym[j] += commvec[i].nEverSymptomatic[j];
      }
    }
    for (int j=0; j<TAG; j++)
      out << "," << nsym[j];
    for (int j=0; j<TAG; j++)
      out << "," << ncsym[j];
    for (int j=0; j<TAG; j++)
      out << "," << ninf[j];
    for (int j=0; j<TAG; j++)
      out << "," << ncinf[j];
    out << endl;
  }
}

/*
 * Create final report in the file `Summary'
 */
void EpiModel::summary(void) {
  ostringstream oss;
  ofstream &outfile = *sumfile;

  {
    outfile << "Corvid version: " << nVersionMajor << "." << nVersionMinor << endl;
    outfile << "Label: " << szLabel << endl;
    outfile << "Population: " << szBaseName << endl;
    outfile << "beta: " << beta << endl;
    outfile << "R0: " << R0 << endl;
    outfile << "Random number generator seed: " << seeddisp << endl;
    outfile << "Tracts: " << nNumTractsTotal << endl;
  }
  outfile << "Communities: " << nNumCommunities << endl;
  outfile << "Families: " << nNumFamilies << endl;
  outfile << "People: " << nNumPerson << endl;
  {
    outfile << "Essential worker fraction: " << fAdultEssentialFraction << endl;
    outfile << "Travel allowed: " << bTravel << endl;
    outfile << "School opening days: ";
    for (int i=0; i<56; i++)
      outfile << nSchoolOpeningDays[i] << ",";
    outfile << endl;

    if (nSeedInfectedTractFIPS>0 || nSeedInfectedCountyFIPS>0 || nSeedInfectedStateFIPS>0)
      outfile << "Seeded " << nSeedInfectedNumber << " infected people in tract " << nSeedInfectedStateFIPS << "," << nSeedInfectedCountyFIPS << "," << nSeedInfectedTractFIPS << endl;
    else
      outfile << "Seeded " << nSeedInfectedNumber << " infected people" << endl;
    outfile << "Seeded (n/10000) at airports: " << nSeedAirports << endl;
    if (bSeedDaily)
      outfile << "Seeded daily" << endl;
    else
      outfile << "Seeded once" << endl;
    outfile << "Protection from pre-existing immunity: " << fPreexistingImmunityProtection << endl;
    outfile << "Pre-existing immunity fraction by age: ";
    for (int i=0; i<TAG; i++)
      outfile << fPreexistingImmunityFraction[i] << ",";
    outfile << endl;
    outfile << "Baseline VES by age: ";
    for (int i=0; i<TAG; i++)
      outfile << fBaselineVESByAge[i] << ",";
    outfile << endl;

    if (ePrevaccinationStrategy==PREVACCINATE)
      outfile << "Prevaccination strategy: prevaccination" << endl;
    else if (ePrevaccinationStrategy==PRIMEBOOSTRANDOM)
      outfile << "Prevaccination strategy: prime/boost random" << endl;
    else if (ePrevaccinationStrategy==PRIMEBOOSTSAME)
      outfile << "Prevaccination strategy: prime/boost same" << endl;
    else
      outfile << "Prevaccination strategy: none" << endl;
    if (eVaccinationStrategy==RINGVACTRACT)
      outfile << "Reactive vaccination strategy: tract" << endl;
    else if (eVaccinationStrategy==RINGVACCOUNTY)
      outfile << "Reactive vaccination strategy: county" << endl;
    else if (eVaccinationStrategy==MASSVAC)
      outfile << "Reactive vaccination strategy: mass" << endl;
    else 
      outfile << "Reactive vaccination strategy: none" << endl;
    outfile << "Vaccine priorities : ";
    for (int i=0; i<PRIORITY_LAST; i++)
      outfile << (int)nVaccinePriorities[i] << ",";
    outfile << endl;
    if (nPriorityChangeTime>0) {
      outfile << "Priority change time: " << nPriorityChangeTime << endl;
      outfile << "Second vaccine priorities: ";
      for (int i=0; i<PRIORITY_LAST; i++)
	outfile << (int)nVaccinePriorities2[i] << ",";
      outfile << endl;
    }

    outfile << "Response threshold: " << fResponseThreshold << endl;
    outfile << "Response delay: " << nTriggerDelay << endl;
    outfile << "Ascertainment delay: " << nAscertainmentDelay << endl;
    outfile << "Ascertainment fraction: " << fSymptomaticAscertainment << endl;
    if (bTrigger)
      outfile << "Reactive strategies deployed on day: " << (int)(nTriggerTime/2) << endl;
    else 
      outfile << "Reactive strategies deployed at time: NA" << endl;
    if (nSchoolClosurePolicy==0)
      outfile << "School closure policy: none" << endl;
    else if (nSchoolClosurePolicy==1)
      outfile << "School closure policy: all" << endl;
    else if (nSchoolClosurePolicy==2)
      outfile << "School closure policy: bytract" << endl;
    else
      outfile << "School closure policy: " << nSchoolClosurePolicy << endl;
    outfile << "School closure days: " << nSchoolClosureDays << endl;
    outfile << "Voluntary isolation: " << fVoluntaryIsolationCompliance << endl;
    outfile << "Ascertained isolation: " << fAscertainedIsolationCompliance << endl;
    outfile << "Quarantine: " << fQuarantineCompliance << endl;
    outfile << "Quarantine duration: " << nQuarantineLength << endl;
    outfile << "Liberal leave compliance: " << fLiberalLeaveCompliance << endl;
    outfile << "Liberal leave duration: " << nLiberalLeaveDuration << endl;
    outfile << "Work from home compliance: " << fWorkFromHomeCompliance << endl;
    outfile << "Work from home duration: " << nWorkFromHomeDuration << endl;
    outfile << "Community contact reduction: " << fCommunityContactReduction << endl;
    outfile << "Community contact reduction duration: " << nCommunityContactReductionDuration << endl;
    outfile << "Antiviral policy: ";
    if (eAVPolicy==NOAV)
      outfile << "none" << endl;
    else if (eAVPolicy==HHTAP)
      outfile << "HHTAP" << endl;
    else if (eAVPolicy==HHTAP100)
      outfile << "HHTAP100" << endl;
    else if (eAVPolicy==FULLTAP)
      outfile << "FULLTAP" << endl;
    else if (eAVPolicy==TREATMENTONLY)
      outfile << "treatmentonly" << endl;
    else
      outfile << "invalid strategy" << endl;
    if (nNumTAPDone>0)
      outfile << "Number of cases TAP: " << nNumTAPDone << endl;
    outfile << "Vaccination fraction: " << fVaccinationFraction << endl;
    outfile << "Vaccines initially available: ";
    if (nNumVaccineTypes>0)
      for (int i=0; i<nNumVaccineTypes; i++)
	outfile << nVaccineInitialSupply[i] << ",";
    else
	outfile << "0";
    outfile << endl;
    outfile << "Vaccines deliverable per day: " << nVaccineDailyLimit << endl;
    outfile << "Antivirals available: " << nAVTotalLimit << endl;
    outfile << "Antivirals deliverable per day: " << nAVDailyLimit << endl;
  
    for (int i=0; i<nNumVaccineTypes; i++) {
      unsigned int sum=0;
      for (int day=0; day<MAXRUNLENGTH; day++)
	sum += vaccineproductionschedule[i][day];
      if (sum>0) {
	outfile << "vaccine " << i << " production/day: ";
	for (int day=0; day<MAXRUNLENGTH; day++)
	  outfile << vaccineproductionschedule[i][day] << ",";
	outfile << endl;
      }
    }
  }

  {
    outfile << "Vaccines used: ";
    if (nNumVaccineTypes>0)
      for (int i=0; i<nNumVaccineTypes; i++)
	outfile << nNumVaccineDosesUsed[i] << ",";
    else
      outfile << "0";
    outfile << endl;
    outfile << "Antivirals used: " << nNumAntiviralsUsed << endl;
    for (int i=0; i<nNumVaccineTypes; i++) {
      outfile << "VEs " << i << ": " << VaccineData[i].VEs << endl;
      outfile << "VEi " << i << ": " << VaccineData[i].VEi << endl;
      outfile << "VEp " << i << ": " << VaccineData[i].VEp << endl;
      outfile << "vaccine " << i << " doses/person: " << VaccineData[i].nNumDoses << endl;
      if (VaccineData[i].nNumDoses>1)
	outfile << "vaccine " << i << " boost on day: " << VaccineData[i].nBoostDay << endl;
      outfile << "vaccine " << i << " restrictions (infants, pre-school age children, school age children, young adults, older adults, elderly, pregnant): " << VaccineData[i].fInfants << "," << VaccineData[i].fAge[0] << "," << VaccineData[i].fAge[1] << "," << VaccineData[i].fAge[2] << "," << VaccineData[i].fAge[3] << "," << VaccineData[i].fAge[4] << "," << VaccineData[i].bPregnant << endl;
      outfile << "vaccine " << i << " efficacy buildup:";
      for (int j=0; j<=VACCEFFLENGTH; j++)
	outfile << " " << VaccineData[i].vacceff[j];
      outfile << endl;
    }

    outfile << "Relative vaccine efficacy by age: ";
    for (int i=0; i<TAG; i++)
      outfile << fVaccineEfficacyByAge[i] << ",";
    outfile << endl;
    outfile << "Needs boost by age and school: ";
    for (int i=0; i<TAG+3; i++)
      outfile << bVaccineBoostByAge[i] << ",";
    outfile << endl;
    outfile << "AVEs: " << AVEs << endl;
    outfile << "AVEp: " << AVEp << endl;
    outfile << "AVEi: " << AVEi << endl;

    for (int i=0; i<nNumVaccineTypes; i++)
      if (nNumVaccineDosesUsedHistory[i].back()>0) {
	outfile << "Cumulative vaccines used (daily) " << i << ": ";
	for (vector< unsigned int >::iterator it = nNumVaccineDosesUsedHistory[i].begin();
	     it != nNumVaccineDosesUsedHistory[i].end();
	     it++)
	  outfile << *it << ",";
	outfile << endl;
      }
    if (nNumAntiviralsUsedHistory.size()>0 && nNumAntiviralsUsed>0) {
      outfile << "Cumulative antivirals used (daily): ";
      for (vector< unsigned int >::iterator it = nNumAntiviralsUsedHistory.begin();
	   it != nNumAntiviralsUsedHistory.end();
	   it++)
	outfile << *it << ",";
      outfile << endl;
    }
    outfile << "Number symptomatic (daily): ";
    for (vector< unsigned int >::iterator it = nNumSymptomaticHistory.begin();
	 it != nNumSymptomaticHistory.end();
	 it++)
      outfile << *it << ",";
    outfile << endl;
    outfile << "Cumulative symptomatic (daily): ";
    for (vector< unsigned int >::iterator it = nNumCumulativeSymptomaticHistory.begin();
	 it != nNumCumulativeSymptomaticHistory.end();
	 it++)
      outfile << *it << ",";
    outfile << endl;
  }

  // count number of symptomatic people
  unsigned int nTotalSym[TAG];	// number of symptomatic by age group
  unsigned int nTotalInf[TAG];	// number of infected by age group
  unsigned int nTotal[TAG];	// number of people
  memset(nTotalSym, 0, sizeof(unsigned int)*TAG);
  memset(nTotalInf, 0, sizeof(unsigned int)*TAG);
  memset(nTotal, 0, sizeof(unsigned int)*TAG);
  for (vector< Tract >::iterator it = tractvec.begin();
       it != tractvec.end();
       it++) {
    Tract &t = *it;
    int nsym[TAG];	// number of infected by age group
    int ninf[TAG];	// number of infected by age group
    int ntot[TAG];	// number of people
    memset(nsym, 0, sizeof(int)*TAG);
    memset(ninf, 0, sizeof(int)*TAG);
    memset(ntot, 0, sizeof(int)*TAG);
    int totalpop = 0;
    for (unsigned int i=t.nFirstCommunity; i<t.nLastCommunity; i++) {
      totalpop += commvec[i].nNumResidents;
      for (int j=0; j<TAG; j++) {
	nsym[j] += commvec[i].nEverSymptomatic[j];
	ninf[j] += commvec[i].nEverInfected[j];
	ntot[j] += commvec[i].nNumAge[j];
      }
    }
    int ni=0, ns=0, nt=0;
    for (int j=0; j<TAG; j++) {
      ns += nsym[j];
      ni += ninf[j];
      nt += ntot[j];
      nTotalSym[j]+=nsym[j];
      nTotalInf[j]+=ninf[j];
      nTotal[j]+=ntot[j];
    }
  }

  unsigned int nVacSymptomatic[TAG*NUMVACCINES];	// number of vaccinated symptomatic people
  unsigned int nUnvacSymptomatic[TAG];	// number of unvaccinated symptomatic people
  unsigned int nVac[TAG*NUMVACCINES];	// number of vaccinated people
  unsigned int nUnvac[TAG];	// number of unvaccinated people
  unsigned int nHighRiskSymptomatic[TAG];	// number of high risk symptomatic people
  unsigned int nPregnantSymptomatic[TAG];	// number of pregnant symptomatic people
  unsigned int nHighRiskPregnantSymptomatic[TAG];	// number of high risk and pregnant symptomatic people
  for (int j=0; j<TAG; j++)
    nUnvacSymptomatic[j] = nUnvac[j] = nHighRiskSymptomatic[j] = nPregnantSymptomatic[j] = nHighRiskPregnantSymptomatic[j] = 0;
  for (int j=0; j<TAG*NUMVACCINES; j++)
    nVac[j] = nVacSymptomatic[j] = 0;
  for (vector< Person >::iterator it = pvec.begin(); 
       it != pvec.end();
       it++) {
    Person &p = *it;
    if (!isSusceptible(p) && getWillBeSymptomatic(p) && 
	(!isInfected(p) ||    // was sick
	 isSymptomatic(p))) { // is sick (do not count people who are not sick yet)
      if (isVaccinated(p))
	nVacSymptomatic[whichVaccine(p)*TAG+p.age]++;
      else
	nUnvacSymptomatic[p.age]++;
      if (isPregnant(p))
	nPregnantSymptomatic[p.age]++;
      if (isHighRisk(p)) {
	nHighRiskSymptomatic[p.age]++;
	if (isPregnant(p))
	  nHighRiskPregnantSymptomatic[p.age]++;
      }
    }
    if (isVaccinated(p))
      nVac[whichVaccine(p)*TAG+p.age]++;
    else
      nUnvac[p.age]++;
  }
  {
    outfile << "Total infection attack rates by age: ";
    for (int j=0; j<TAG; j++)
      outfile << (nTotalInf[j]/(double)nTotal[j]) << ",";
    outfile << endl;
    unsigned int ni=0, ns=0, nt=0;
    for (int j=0; j<TAG; j++) {
      ni += nTotalInf[j];
      ns += nTotalSym[j];
      nt += nTotal[j];
    }
    outfile << "Total infection attack rate: " << (ni/(double)nt) << endl;
    outfile << "Total infected individuals by age: ";
    for (int j=0; j<TAG; j++)
      outfile << nTotalInf[j] << ",";
    outfile << endl;

    outfile << "Total symptomatic attack rates by age: ";
    for (int j=0; j<TAG; j++)
      outfile << (nTotalSym[j]/(double)nTotal[j]) << ",";
    outfile << endl;
    outfile << "Total symptomatic attack rate: " << (ns/(double)nt) << endl;
    outfile << "Total symptomatic individuals by age: ";
    for (int j=0; j<TAG; j++)
      outfile << nTotalSym[j] << ",";
    outfile << endl;

    outfile << "Total symptomatic unvaccinated individuals by age: ";
    for (int j=0; j<TAG; j++) {
      outfile << nUnvacSymptomatic[j] << ",";
    }
    outfile << endl;
    for (int i=0; i<nNumVaccineTypes; i++)
      if (nVac[i*TAG] || nVac[i*TAG+1] || nVac[i*TAG+2] || nVac[i*TAG+3] || nVac[i*TAG+4]) {
	outfile << "Total symptomatic vaccinated individuals by age (vaccine " << i << "): ";
	for (int j=0; j<TAG; j++)
	  outfile << nVacSymptomatic[i*TAG+j] << ",";
	outfile << endl;
      }
    outfile << "Total unvaccinated individuals by age: ";
    for (int j=0; j<TAG; j++)
      outfile << nUnvac[j] << ",";
    outfile << endl;
    for (int i=0; i<nNumVaccineTypes; i++)
      if (nVac[i*TAG] || nVac[i*TAG+1] || nVac[i*TAG+2] || nVac[i*TAG+3] || nVac[i*TAG+4]) {
	outfile << "Total vaccinated individuals by age (vaccine " << i << "): ";
	for (int j=0; j<TAG; j++)
	  outfile << nVac[i*TAG+j] << ",";
	outfile << endl;
      }
    if (fHighRiskFraction[0]>0.0 || fHighRiskFraction[1]>0.0 || fHighRiskFraction[2]>0.0 || fHighRiskFraction[3]>0.0 || fHighRiskFraction[4]>0.0) {
      outfile << "High risk fraction by age: ";
      for (int j=0; j<TAG; j++)
	outfile << fHighRiskFraction[j] << ",";
      outfile << endl;
      outfile << "Symptomatic high risk individuals by age: ";
      for (int j=0; j<TAG; j++)
	outfile << nHighRiskSymptomatic[j] << ",";
      outfile << endl;
    }
    if (fPregnantFraction[0]>0.0 || fPregnantFraction[1]>0.0 || fPregnantFraction[2]>0.0 || fPregnantFraction[3]>0.0 || fPregnantFraction[4]>0.0) {
      outfile << "Pregnant fraction by age: ";
	for (int j=0; j<TAG; j++)
	  outfile << fPregnantFraction[j] << ",";
      outfile << endl;
      outfile << "Symptomatic pregnant individuals by age: ";
      for (int j=0; j<TAG; j++)
	outfile << nPregnantSymptomatic[j] << ",";
      outfile << endl;
      if (fHighRiskFraction[0]>0.0 || fHighRiskFraction[1]>0.0 || fHighRiskFraction[2]>0.0 || fHighRiskFraction[3]>0.0 || fHighRiskFraction[4]>0.0) {
	outfile << "Symptomatic high risk and pregnant individuals by age: ";
	for (int j=0; j<TAG; j++)
	  outfile << nHighRiskPregnantSymptomatic[j] << ",";
	outfile << endl;
      }
    }
    outfile << "Total individuals by age: ";
    for (int j=0; j<TAG; j++)
      outfile << nTotal[j] << ",";
    outfile << endl;
  }

    outfile.close(); 
}

void EpiModel::outputIndividuals(void) {
  if (bIndividualsFile) {
    ostream &out = *individualsfile;
    out << "id,age,familyid,homecomm,homeneighborhood,daycomm,dayneighborhood,workplace,infectedtime,incubationdays,sourceid,sourcetype,vacstatus" << endl;
    for (vector< Person >::iterator it = pvec.begin();
       it != pvec.end();
       it++) {
      Person &p = *it;
      out 
	  << p.id << "," << (int)p.age << "," << p.family << ","
	  << p.nHomeComm << "," << (int)p.nHomeNeighborhood << "," 
	  << p.nDayComm <<  "," << (int)p.nDayNeighborhood << "," << (int)p.nWorkplace << "," 
	  << p.nInfectedTime << "," << p.nIncubationDays << "," << p.sourceid << "," << (int) p.sourcetype << "," <<  (int) isVaccinated(p) << endl; 
    }
    (*individualsfile).close();
  }
}

/*
 * Randomly infect individuals?
 */
void EpiModel::seedinfected(void) {
  // seeding infected persons
  if (nSeedInfectedNumber>0) {
    if (nSeedInfectedTractFIPS>0 || nSeedInfectedCountyFIPS>0 || nSeedInfectedStateFIPS>0) {
      for (vector< Tract >::iterator it = tractvec.begin();
	   it != tractvec.end();
	   it++) {
	Tract &t = *it;
	if (t.fips_state==nSeedInfectedStateFIPS &&
	    t.fips_county==nSeedInfectedCountyFIPS &&
	    t.fips_tract==nSeedInfectedTractFIPS) {
	  for (int i=0; i<nSeedInfectedNumber; i++) { // infect people in this tract
	    long c = t.nFirstCommunity;
	    if (t.nLastCommunity-t.nFirstCommunity>1)
	      c += get_rand_uint32 % (t.nLastCommunity-t.nFirstCommunity);
	    long pid = commvec[c].nFirstPerson + get_rand_uint32 % commvec[c].nNumResidents;
	    infect(pvec[pid]);
	    pvec[pid].sourcetype = FROMSEED;
	  }
	}
      }
    } else {
      // seed infected people across whole population
      for (int i=0; i<nSeedInfectedNumber; i++) { // infect people
	long pid = get_rand_uint32 % nNumPerson;
	if (isSusceptible(pvec[pid])) {
	  infect(pvec[pid]);
	  pvec[pid].sourcetype = FROMSEED;
	}
      }
    }
  }
  // airport seeding is independent of other seeding
  if (nSeedAirports>0) {
    for (vector< Tract >::iterator it = tractvec.begin();
	 it != tractvec.end();
	 it++) {
      Tract &t = *it;
      for (int hub=0; hub<NUMBER_HUBS; hub++) {
	if (t.fips_state==FIPS_hubs[hub]/1000 &&
	    t.fips_county==FIPS_hubs[hub]%1000) {

	  for (; it != tractvec.end() && (*it).fips_state==t.fips_state && (*it).fips_county==t.fips_county; it++)
	    ;
	  --it;
	  unsigned int nFirst = commvec[t.nFirstCommunity].nFirstPerson;
	  unsigned int nRange = commvec[(*it).nLastCommunity-1].nLastPerson-nFirst;
	  
	  int num = bnldev(nSeedAirports/(float)10000, Size_hubs[hub]/365);
	  for (int i=num; i>0; i--) { // infect people in this county
	    unsigned long pid = nFirst + get_rand_uint32 % nRange;
	    assert(pid<pvec.size());
	    infect(pvec[pid]);
	    pvec[pid].sourcetype = FROMAIRPORT;
	  }
	}
      }
    }
  }
}

/*
 * Actions to take just before the model run
 */
void EpiModel::prerun(void) {
  // open schools as appropriate
  for (vector< Tract >::iterator it = tractvec.begin();
       it != tractvec.end();
       it++) {
    Tract &t = *it;
    if (nSchoolOpeningDays[t.fips_state-1]>0) {
      for (int i=0; i<9; i++)
	setSchoolClosed(t,i);// school not open yet in this state
    }
  }

  // pre-vaccination and priming
  if (fVaccinationFraction>0.0 &&
      (ePrevaccinationStrategy==PREVACCINATE ||
       ePrevaccinationStrategy==PRIMEBOOSTRANDOM || 
       ePrevaccinationStrategy==PRIMEBOOSTSAME)) {
    // pre-vaccinate or pre-prime those who want it
    vector< Person >::iterator pend=pvec.end();
    for (vector< Person >::iterator it = pvec.begin(); 
	 it != pend;
	 it++) {
      Person &p = *it;
      if (p.nVaccinePriority>0 &&
	  (ePrevaccinationStrategy!=PRIMEBOOSTRANDOM ||
	   (ePrevaccinationStrategy==PRIMEBOOSTRANDOM && get_rand_double<fVaccinationFraction)))
	p.bWantVac=true;
    }
    for (int vacnum=0; vacnum<nNumVaccineTypes; vacnum++)
      if (VaccineData[vacnum].nNumDoses>0) {
	unsigned int nDoses = (ePrevaccinationStrategy==PREVACCINATE?VaccineData[vacnum].nNumDoses:1); // for prevaccination, give all doses; For priming, just give 1.
	int nDay = (ePrevaccinationStrategy==PREVACCINATE?VACCEFFLENGTH:VaccineData[vacnum].nBoostDay); // how many days ago to pre-vaccinate or pre-prime?
	double vFrac = 0.0;
	// demand count for THIS vaccine
	unsigned int wantthisvac[PRIORITY_LAST];
	memset(wantthisvac, 0, PRIORITY_LAST*sizeof(unsigned int));
	for (vector< Person >::iterator it = pvec.begin(); 
	     it != pend;
	     it++) {
	  Person &p = *it;
	  if (p.bWantVac && p.bVaccineEligible[vacnum] && p.nVaccinePriority>0)
	    wantthisvac[p.nVaccinePriority]++;
	}
	for (int priority=1; priority<=nHighestPriority; priority++) {
	  if (wantthisvac[priority]>0) {
	    vFrac = nVaccineSupply[vacnum]/(double)wantthisvac[priority];
	    for (vector< Person >::iterator it = pvec.begin(); 
		 it != pend;
		 it++) {
	      Person &p = *it;
	      if (p.bWantVac && p.bVaccineEligible[vacnum] && p.nVaccinePriority>0) {
		if (vFrac>=1.0 || get_rand_double < vFrac) {
		  p.bWantVac = false;
		  if (nVaccineSupply[vacnum]>=nDoses)
		    nVaccineSupply[vacnum]-=nDoses;
		  else
		    nVaccineSupply[vacnum]=0;
		  nNumVaccineDosesUsed[vacnum]+=nDoses;
		  setWhichVaccine(p, vacnum);
		  setVaccinated(p);
		  if (needsBoost(p)) {
		    if (VaccineData[vacnum].nNumDoses>1) {
		      clearNeedsBoost(p);
		    } else if (ePrevaccinationStrategy==PREVACCINATE && nVaccineSupply[vacnum]>0) {
		      nVaccineSupply[vacnum]--;
		      nNumVaccineDosesUsed[vacnum]++;
		    }
		  }
		  if (ePrevaccinationStrategy==PREVACCINATE)
		    setBoosted(p);
		  p.vday = nDay;
		}
	      }
	    }
	  }
	}
      }

    // clear bWantVac status of all people
    for (vector< Person >::iterator it = pvec.begin(); 
	 it != pvec.end();
	 it++) {
      Person &p = *it;
      p.bWantVac = false;
      if ((ePrevaccinationStrategy==PRIMEBOOSTSAME && !isVaccinated(p)) || // only vaccinated people will get a boost
	  (ePrevaccinationStrategy==PRIMEBOOSTRANDOM && get_rand_double>=fVaccinationFraction)) // random people will get a dose
	p.nVaccinePriority = 0;
    }
  }
  if (!bSeedDaily)
    seedinfected();
}

/* 
 * run - runs the simulation for the desired number of days
 *   nTimer: initialized to 0 in the constructor
 *   nRunLength: total number of days requested
 */
void EpiModel::run(void) {
  prerun();
  while(nTimer<nRunLength*2) {
    if (nLogFileInterval>0 && (int)(nTimer/2)%nLogFileInterval==0)
      log();
    if (bSeedDaily) {
      seedinfected();
    }
    day();
    if (bTravel)
      travel_end();
    nTimer++;
    night();
    if (bTravel)
      travel_start();

    response();
    nTimer++; 

    unsigned int nNumSymptomatic=0;
    unsigned int nNumCumulativeSymptomatic=0;
    for (vector< Community >::iterator it = commvec.begin();
	 it != commvec.end();
	 it++) {
      for (int j=0; j<TAG; j++) {
	nNumCumulativeSymptomatic += (*it).nEverSymptomatic[j];
	nNumSymptomatic += (*it).nsym[j];
      }
    }
    nNumCumulativeSymptomaticHistory.push_back(nNumCumulativeSymptomatic);
    nNumSymptomaticHistory.push_back(nNumSymptomatic);
  }
  if (logfile)
    (*logfile).close();
  summary();
  outputIndividuals();
}
