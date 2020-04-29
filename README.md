
# Purpose

This is the source code of the hybrid genetic search (HGS) for the TDCARP, described in 
"Vidal, T., Martinelli, R., Pham, T. A., & Hà, M. H. (2020). Arc routing with time-dependent travel times and paths"

For convenience, the calculation of the quickest path is done in a separate code. The metaheuristic reads the result of the quickest path calculation as a TDSP file, and finds near-optimal itineraries to service the requests.
This implementation is an extension of the HGS for capacitated arc routing problems available at "https://github.com/vidalt/HGS-CARP".

# Running

```
Usage:
   ./gencarp tdsp_file_path [options]
Available options:
  -t             CPU time in seconds (defaults to 300s).
  -sol           File where to output the solution statistics (defaults to the instance file name prepended with 'sol-').
  -seed          Seed for initialization of the RNG
```

Example: `./gencarp ../Instances/TDCARP-withSP/Type_M/C01.tdsp -sol solution.txt -t 300 -seed 1`

