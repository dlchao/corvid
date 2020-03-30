#!/bin/bash
#script for running total school closures (non-reactive)
#timing and duration of closure are swept
for r0 in 2.6; do
    for responseday in 60 74 100; do
	for duration in 14 42 180; do
	    printr0=${r0//./}
	    cp template-config-seattle-totalschoolclosure configschool
	    sed -i "s/PARAMR0/$r0/g" configschool
	    sed -i "s/PRINTR0/$printr0/g" configschool
	    sed -i "s/DURATION/$duration/g" configschool
	    sed -i "s/RESPONSEDAY/$responseday/g" configschool
	    ./corvid configschool
	done
    done
done

exit 0
