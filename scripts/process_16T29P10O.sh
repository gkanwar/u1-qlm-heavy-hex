#!/bin/bash

shape=16T29P10O
for Kp in 0.40 0.70 2.00; do
    python scripts/analyze_specific_geom.py --KP ${Kp} --geom ${shape} --E0=26
done
