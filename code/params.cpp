/*
 * params.cpp
 * constants for Corvid
 */

#include <climits>
#include "params.h"

const double fAG0InfantFraction=0.09923; // fraction of those in age group 0 who are <6 months

// contact probabilities
const double cpw= 0.05; // workplace
const double cps[10] = {0.0, 0.0252, 0.03, 0.0348, 0.0348, 0.12, 0.12, 0.12, 0.12, 0.28}; // school, 10 values for:
//      none, high school, middle school, elementary (2), day care (4), play group
const double cpcm[TAG] = {0.0000109, 0.0000326, 0.000087, 0.000087, 0.000174}; // community
const double cpnh[TAG] = {0.0000435, 0.0001305, 0.000348, 0.000348, 0.000696}; // neighborhood
const double cpfc[TAG] = {0.8, 0.8, 0.37, 0.37, 0.37};         // family from children
const double cpfa[TAG] = {0.25, 0.25, 0.37, 0.37, 0.37};       // family from adults
const double cphcc[TAG] = {0.08, 0.08, 0.037, 0.037, 0.037};   // household cluster from children
const double cphca[TAG] = {0.025, 0.025, 0.037, 0.037, 0.037}; // household cluster from adults

// withdraw probabilities
const double withdrawprob[3][WITHDRAWDAYS] = {
  {0.304, 0.575, 0.324},  // preschool children
  {0.203, 0.498, 0.375},  // school children
  {0.100, 0.333, 0.167}	  // adults
};

// self isolation probabilities by age
const double isolationprob[TAG] = {0.8,0.75,0.5,0.5,0.5};

const int nAntiviralCourseSize = 10; // number of pills in one antiviral course (1 tablet/day for prophylaxis and 2 tablets/day for treatment). based on oseltamivir.
const double fStopAntiviralTwoPills = 0.0; // probability that individuals taking antivirals stop after exactly two pills (parameter was used for Oseltamivir)

// symptomatic probability
double fBaseSymptomaticProb=0.5; // half of infections are symptomatic

// ratio of symptomatic:asymptomatic infectiousness
double fRelativeSymptomaticInfectiousness = 2.0;

// cdf of duration for incubation period in days
//const double incubationcdf[14] = {0.0128885,0.06283671,0.1533217,0.2772347,0.4202624,0.5651,0.6961202,0.8029434,
//				  0.8817976,0.9346127,0.9667359,0.9844785,0.9933733,1.0}; // Weibull from Hellewell et al 2020
const double incubationcdf[INCUBATIONMAX] = {0.000053,0.013218,0.105699,0.287229,0.488966,0.658552,0.781508,0.863625,
				  0.91598,0.948511,0.968458,0.980622,0.988035,1.0}; // log normal from Lauer et al 2020

// viral load trajectories (infectiousness)
const double basevload[VLOADNSUB][VLOADNDAY] =
  {{0.00015087,0.00394974,0.02116544,0.05230587,0.08463726,0.1065278,0.114426,
0.110691,0.09968273,0.08537833,0.07055667,0.05682167,0.0449103,0.03501563,
    0.02703378,0.02072588,0.01581306,0.01202629,0.0091288,0.00692299}, // log normal
   {0.00015087,0.00394974,0.02116544,0.05230587,0.08463726,0.1065278,0.114426,
0.110691,0.09968273,0.08537833,0.07055667,0.05682167,0.0449103,0.03501563,
    0.02703378,0.02072588,0.01581306,0.01202629,0.0091288,0.00692299}, //  just a copy of trajectory 1
   {0.00015087,0.00394974,0.02116544,0.05230587,0.08463726,0.1065278,0.114426,
0.110691,0.09968273,0.08537833,0.07055667,0.05682167,0.0449103,0.03501563,
    0.02703378,0.02072588,0.01581306,0.01202629,0.0091288,0.00692299}}; //  just a copy of trajectory 2

// travel data
double const travel_pr[TAG] = { .0023, .0023, .0050, .0053, .0028 }; // age-specific travel probability per day
double const travel_length_cdf[12] = { 0.239, 0.406, 0.574, 0.741,
				       0.787, 0.834, 0.880, 0.926,
				       0.945, 0.963, 0.982, 1.0 }; // cdf of length of trip in days (nights away from home)

const unsigned int FIPS_hubs[] = {36061, // JFK
				  6037,  // LAX
				  12086, // MIA
				  17031, // ORD
				  34017, // EWR
				  13121, // ATL
				  6075,  // SFO
				  48201, // IAH
				  51059, // IAD
				  48113, // DFW
				  26163, // DTW
				  42045, // PHL
				  25025, // BOS
				  12011, // FLL
				  53033};// SEA
const unsigned int Size_hubs[] = {21842544, 17019166, 15509279, 11375367, 
				  10812993, 9166055,  8648219,  7627942,
				  5893142,  4872207,  3887481,  3734127,
				  3673748,  3062384,  2766576};
// from US International Air Passenger and Freight Statistics, June 2008
