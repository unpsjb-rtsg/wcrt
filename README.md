# WCRT
Evaluation of exact schedulability methods.

## Overview
The following exact schedulability methods are implemented:
+ RTA, in [Improved Response-Time Analysis Calculations](http://doi.ieeecomputersociety.org/10.1109/REAL.1998.739773).
+ RTA2, in [Reduced computational cost in the calculation of worst case response time for real time systems](http://sedici.unlp.edu.ar/handle/10915/9654).
+ RTA3, in [Computational Cost Reduction for Real-Time Schedulability Tests Algorithms](http://ieeexplore.ieee.org/document/7404899/).
+ HET, in [Schedulability Analysis of Periodic Fixed Priority Systems](http://ieeexplore.ieee.org/document/1336766/).

Implementations in C and Python are provided.

## Python version
Requires Python3. This script outputs its results into a HDFS store.

### Dependencies
+ Pandas 

## C version
To compile `wcrt-test-sim.c`:

     # gcc -o wcrt-test-sim wcrt-test-sim.c -Wall -I/usr/include/libxml2 -Wall -lgsl -lxml2

### Dependencies
+ libxml2
+ gsl - GNU Scientific Library

# License
This software is licensed under the MIT License. A copy of the license can be found in the `LICENSE` file.
