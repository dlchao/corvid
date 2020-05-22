#!/bin/bash
# Figure 7A-B
# Schools close on day 30, re-open on day 180

#ll/wfh
./corvid config-seattle26-layered1-30-none-150-365-0.5-0-0-0-none-0-14-1
#ll/wfh+school 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0-0-0-none-0-14-1
#ll/wfh+school+shelter 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0-0-none-0-14-1
#ll/wfh+school+shelter+home isolate 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-0-none-0-14-1
#ll/wfh+school+shelter+self isolate 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-0-treatmentonly-0-14-1
#ll/wfh+school+shelter+home isolate+quarantine 
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-0-none-0.8-14-1
./corvid config-seattle26-layered1-30-all-150-365-0.5-0.3-0.5-0-none-0.8-14-1
