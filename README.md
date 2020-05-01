
# Description

This is the source code of the hybrid genetic search (HGS) for the TDCARP, described in 
"Vidal, T., Martinelli, R., Pham, T. A., & Ha, M. H. (2020). Arc routing with time-dependent travel times and paths"

The code permits to perform the calculations from a known network. At this point, it does not include data acquisition routines, geocoding routines, and interfaces.
It is entirely in C++. For convenience, the calculation of the continuous quickest path is done in a separate folder named `Program-SP` whereas the metaheuristic is located in the `Program` folder.
The metaheuristic is an extension of the HGS for capacitated arc routing problems available at "https://github.com/vidalt/HGS-CARP", with new route evaluation operations and lower-bounds for move evaluations.

# Using the Metaheuristic

* Building
```
cd Program
make
```

* Running:
```
Usage:
   ./gencarp file_path.tdsp [options]
Available options:
  -t             CPU time in seconds (defaults to 300s).
  -sol           File where to output the solution statistics (defaults to the instance file name prepended with 'sol-').
  -seed          Seed for initialization of the RNG
```

Example: `./gencarp ../Instances/TDCARP-withSP/Type_M/C01.tdsp -sol solution.txt -t 300 -seed 1`

# Using the Quickest Path algorithm

* Building
```
cd Program-SP
make
```

* Running:
```
Usage:
./TDSPGenerator file_path.dat
```
Example: `./TDSPGenerator ../Instances/TDCARP-withoutSP/Type_H/C01.dat`
This command creates a new tdsp file named `./C01.tdsp`
