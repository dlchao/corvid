/* main program for Corvid
 *
 * Dennis Chao
 * May 2020
 */

#include "epimodel.h"
#include "epimodelparameters.h"

int main(int argc, char *argv[]) {
  char *configname=NULL;
  if (argc==2)
    configname = argv[1];
  else {
    cerr << "Usage: " << argv[0] << " configfilename" << endl;
    exit(-1);
  }
  EpiModelParameters parms(configname);
  EpiModel model(parms);
  model.run();
  return 0;
}
