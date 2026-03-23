# EcoEvolutionSim: Architecture and Implementation

This document provides a detailed technical breakdown of the simulation's design and logic following the structural refactor.

---

## 1. Data Architecture: Structure of Arrays (SoA)

Unlike many agent-based models that store a list of agent objects, this simulation uses a **Structure of Arrays (SoA)** pattern.

*   **Implementation:** All agent properties (positions, velocities, traits, energy) are stored in individual `Vector{Float32}` blocks within the `Agents` struct.
*   **Performance:** This layout is highly optimized for modern CPUs. It ensures **cache locality**, meaning that when the CPU processes a list of $x$ coordinates, it can pre-load the next coordinates into its fast cache. It also enables **SIMD** (Single Instruction, Multiple Data) optimizations.
*   **Memory Efficiency:** The use of `Float32` instead of `Float64` reduces the memory footprint by 50% and speeds up arithmetic on most hardware.

---

## 2. Spatial Engine: Parallel Cell Grid

To maintain high performance as the population grows, the world is partitioned into a **Spatial Grid**.

1.  **Binning:** Every agent is assigned to a grid cell based on its coordinates.
2.  **Prefix Sum Grid:** During each step, the `build_cell_grid!` function rebuilds the spatial index. It uses `Threads.Atomic` operations to ensure that multiple CPU cores can safely count and bin agents simultaneously.
3.  **Local Search:** When calculating interactions, an agent only checks other agents in its own cell and the 8 immediate neighbors. This reduces the complexity from $O(N^2)$ to $O(N)$, allowing the simulation to handle hundreds of thousands of agents.

---

## 3. Modular Registry System

The core of the new architecture is the **Registry System** (`src/registry/`), which decouples the simulation engine from the biological logic.

*   **Trait Registry:** Manages different trait types (e.g., `continuous`, `bounded`).
*   **Fitness Registry:** Manages ecological interaction "plugins."

### How Plugins Work:
1.  The `Config` (from `config.toml`) specifies an interaction type string (e.g., `"default_interaction"`).
2.  During initialization, the engine looks up this string in the `FITNESS_REGISTRY`.
3.  The registry returns a **Constructor** that builds a specific **Kernel Function**.
4.  The engine then stores this function as a pointer (`sim.fitness_fn`) and executes it every step without needing to know the math inside it.

---

## 4. The Simulation Loop (`step!`)

The `step!` function in `src/simulator.jl` is the main execution cycle:

1.  **Spatial Update:** `build_cell_grid!` updates the grid based on the latest positions.
2.  **Ecological Step:** `compute_interactions!` executes the "compiled" `fitness_fn`. This handles competition, mating, and predation.
3.  **Physical Step:** `movement_step!` updates positions and velocities. It supports various strategies like **Langevin** dynamics (simulating momentum and friction) and applies **Reflective Boundaries** to keep agents within the world.

---

## 5. Ecological Logic

The default biology logic is defined in `src/fitness.jl`. It uses Gaussian kernels to determine interaction strength based on trait distance ($d$):

*   **Competition:** $cost = e^{-d^2 / 2\sigma^2}$. Similar traits lead to higher energy loss, simulating niche competition.
*   **Mating:** If trait similarity exceeds a threshold (0.8), agents receive an energy bonus, representing reproductive success or cooperation.
*   **Predation:** If the trait difference $d$ is larger than a threshold, the agent with the higher trait value acts as a predator, gaining energy at the expense of the prey.

---

## 6. Visualization and Observables

The `script/visualize.jl` script uses the **Makie.jl** ecosystem for high-performance rendering.

*   **Observables:** Instead of redrawing the entire frame, the script uses `Observables` (e.g., `points[] = ...`). This allows the GPU to update only the changed data, maintaining high frame rates.
*   **Visual Mapping:** Agents are colored based on their trait values using a colormap (e.g., `:viridis`), and their size is scaled by their current energy level. This allows for real-time observation of emergent population clusters and "species."

---

## 7. Design Benefits

By moving to this modular, registry-based architecture, the simulation becomes a **Scientific Engine**:
*   **Extensibility:** Researchers can add new traits or ecological behaviors by writing a small function and registering it, without ever touching the complex grid or movement code.
*   **Maintainability:** Each component (Movement, Environment, Fitness, Registry) is isolated, making it easier to test and debug.
*   **Performance:** The engine is optimized for multi-core CPUs while keeping the user-facing logic simple and biological.
