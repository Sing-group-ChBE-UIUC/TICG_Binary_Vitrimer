# Brownian Dynamics Simulation with Bond Exchange

This repository contains a C-based hybrid Brownian Dynamics (BD) / Monte Carlo (MC) simulation capable of modeling polymer networks with dynamic bond exchange reactions. 

## Features
* **Theoretically-Informed Coarse Graining:** Polymer segments interact with TICG bonded and nonboneded potentials. 
* **Brownian Dynamics:** Simulates polymer movements using Langevin dynamics.
* **Bond Exchange Reaction:** Models bond exchange reactions using a Monte Carlo Metropolis algorithm.

## Compilation
To compile the simulation, you need a C compiler (like `gcc`). The code relies on the standard math library.

```bash
gcc -o simulation main.c -lm
```

## Authors
Yun-Ju Chen and Charles E. Sing

## Funding Acknowledgements

This work is supported by the United States National Science Foundation (CBET-2029928) and the U.S. Department of Energy, Office of Basic Energy Sciences, Division of Materials Sciences and Engineering, under Award DE-SC0020858. The authors acknowledge the facility and instrumental support from the Materials Research Laboratory and the SCS NMR Laboratory, University of Illinois Urbana-Champaign.
