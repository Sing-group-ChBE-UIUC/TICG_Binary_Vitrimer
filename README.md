# Brownian Dynamics Simulation with Bond Exchange

This repository contains a C-based hybrid Brownian Dynamics (BD) / Monte Carlo (MC) simulation capable of modeling polymer networks with dynamic bond exchange reactions. 

## Features
* **Brownian Dynamics:** Simulates polymer movements using Langevin dynamics.
* **Dynamic Bonding:** Models bond exchange reactions using a Monte Carlo Metropolis algorithm.

## Compilation
To compile the simulation, you need a C compiler (like `gcc`). The code relies on the standard math library.

```bash
gcc -o simulation main.c -lm
