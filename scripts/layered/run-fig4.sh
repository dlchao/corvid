#!/bin/bash
# Figure 7
# test and isolate with delays

#ll/wfh+school+shelter+home isolate 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-5-none-0-14-1
#ll/wfh+school+shelter+self isolate 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-5-treatmentonly-0-14-1
#ll/wfh+school+shelter+home isolate+quarantine 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-5-none-0.8-14-1

#ll/wfh+school+shelter+home isolate 
./corvid config-seattle26-layered1-90-all-90-365-0.5-0.3-0.5-5-none-0-14-1
#ll/wfh+school+shelter+self isolate 
./corvid config-seattle26-layered1-90-all-90-365-0.5-0.3-0.5-5-treatmentonly-0-14-1
#ll/wfh+school+shelter+home isolate+quarantine 
./corvid config-seattle26-layered1-90-all-90-365-0.5-0.3-0.5-5-none-0.8-14-1
