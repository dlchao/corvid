#!/bin/bash
# interventions start on day 30 and never stop

#baseline
./corvid config-seattle26-layered1-30-none-365-365-0-0-0-0-none-0-14-1
#ll/wfh
./corvid config-seattle26-layered1-30-none-365-365-0.5-0-0-0-none-0-14-1
#ll/wfh+school 
./corvid config-seattle26-layered1-30-all-365-365-0.5-0-0-0-none-0-14-1
#ll/wfh+school+shelter 
./corvid config-seattle26-layered1-30-all-365-365-0.5-0.3-0-0-none-0-14-1
#ll/wfh+school+shelter+home isolate 
./corvid config-seattle26-layered1-30-all-365-365-0.5-0.3-0.5-0-none-0-14-1
#ll/wfh+school+shelter+self isolate 
./corvid config-seattle26-layered1-30-all-365-365-0.5-0.3-0.5-0-treatmentonly-0-14-1
#ll/wfh+school+shelter+home isolate+quarantine 
./corvid config-seattle26-layered1-30-all-365-365-0.5-0.3-0.5-0-none-1-14-1
