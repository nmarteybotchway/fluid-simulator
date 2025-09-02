# Fluid Simulation

A real-time 2D fluid simulation in C++ using **SFML**, implementing a stable fluid solver based on Jos Stam's method.  
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

![Fluid Simulation Screenshot](screenshot.png)  <!-- Replace with actual screenshot if available -->

---

## Getting Started

### Prerequisites

- **C++17** or higher
- [SFML](https://www.sfml-dev.org/) (tested with SFML 2.5+)
- [CMake](https://cmake.org/) (for building the project)
- Optional: **OpenMP** support for parallelization

### Building

```bash
git clone <your-repo-url>
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

## Project Structure

### `Fluid.h` / `Fluid.cpp`

Main fluid simulation class implementing:

- `addDensity` / `addVelocity` – Add density or velocity to the fluid at a given cell.
- `advanceSimulation` – Performs a full simulation step: advection, diffusion, and projection.
- `renderDensity` – Updates an SFML texture for visualization.
- Optional: `decayDensity` – Gradually reduces density over time (can be disabled).

### `Grid.h`

Simple 2D float grid wrapper with bounds-checked access, used for density and velocity fields.

### `main.cpp`

Entry point of the application. Handles:

- SFML window creation and rendering loop.
- Mouse input to inject density/velocity.
- Calls `advanceSimulation` and `renderDensity` every frame.