# Fluid Simulation

A real-time 2D fluid simulation in C++ using **SFML**, implementing a stable fluid solver.  
Add density and velocity interactively with your mouse to create swirling fluid effects.

---

## Features

- **Real-time fluid simulation** with advection, diffusion, and projection steps.
- **Interactive input**: Click and drag the mouse to add density and velocity.
- **Smooth rendering** with SFML textures.
- **Parallelized loops** using OpenMP for faster computation.
- **Optional decay** of density for fading effects.

---

## Demo

![Fluid Simulation Screenshot](assets/demo.gif)  <!-- Replace with actual screenshot if available -->

---

## Getting Started

### Prerequisites

- **C++17** or higher
- [SFML](https://www.sfml-dev.org/) (tested with SFML 2.5+)
- [CMake](https://cmake.org/) (for building the project)
- Optional: **OpenMP** support for parallelization

### Building

```bash
git clone https://github.com/nmarteybotchway/fluid-simulator.git
cd fluid_simulation
mkdir build && cd build
cmake ..
cmake --build .
```

### Running

```bash
./fluid_simulation
```

- Left mouse click + drag: add density and velocity.
- Close window to exit

## References

1. **Mike Ash, "Fluid Simulation for Dummies"**  
   https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

2. **OpenMP Documentation – Parallel Programming in C++**  
   [https://www.openmp.org/resources/](https://www.openmp.org/resources/)

3. **SFML (Simple and Fast Multimedia Library)**  
   [https://www.sfml-dev.org/](https://www.sfml-dev.org/)  
   
