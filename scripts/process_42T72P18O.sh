#!/bin/bash

shape=42T72P18O
for Kp in 0.40 0.70 2.00; do
    python scripts/analyze_specific_geom.py --KP ${Kp} --geom ${shape} --E0=68
done
