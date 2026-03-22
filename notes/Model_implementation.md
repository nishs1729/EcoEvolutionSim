# EcoEvolutionSim: Implementation Details

## 1. The Core Model: Phenotype-Space Ecology

The simulation is built on the principle that every individual agent is defined by a **trait vector** (phenotype). Instead of hard-coding specific behaviors or strategies, behaviors emerge from these traits and their interactions with the environment and other agents.

### Current Traits
In the current implementation, we model:
- **Ecological Niche/Resource Preference ($z_1$):** A continuous value in $[0, 1]$ that determines which "resources" an agent is best at using and who it competes with.
- **Energy State:** A dynamic variable representing the "health" or "fitness" of the agent, affected by interactions.

## 2. Spatial Modeling: The Cell Grid

To handle millions of agents efficiently, the world is divided into a **Spatial Grid**.
- **World Size:** Defined in the `Config`.
- **Cell Size:** Determines the resolution of spatial interactions.
- **Complexity:** This reduces the interaction search from $O(N^2)$ to $O(N \times \text{local density})$, which is critical for performance.
- **Implementation:** Uses a "Structure of Arrays" (SoA) for agent data and a prefix-sum-based cell grid for thread-safe parallel processing.

## 3. Movement Strategies

Movement is modeled to be smooth yet diffusive, following various physical and biological patterns.

### RANDOM_WALK
The baseline model where agents move by adding a random displacement in each step.
- **Dynamics:** $x_{t+1} = x_t + \eta$ where $\eta$ is uniform noise.

### LANGEVIN (Smooth Velocity-based)
Simulates particles in a fluid with friction and random kicks.
- **Dynamics:** $v_{t+1} = v_t(1 - \gamma) + \text{noise}$, where $\gamma$ is friction.
- **Visuals:** Produces smooth, gliding motion with momentum.

### CORRELATED_RW (Persistent Motion)
Agents have a "heading" and prefer to move forward, with their heading changing slowly.
- **Dynamics:** $\theta_{t+1} = \theta_t + \text{noise}$; $x_{t+1} = x_t + v \cos(\theta)$.
- **Visuals:** Produces curved paths and "intent-like" movement.

### ACTIVE_BROWNIAN
Similar to the Correlated Random Walk but often used to model self-propelled particles with constant speed.

## 4. Ecological Interactions

Interactions occur when two agents are within the `interaction_radius`. The outcome depends on their trait distance $d = |z_i - z_j|$.

### Competition
Agents with similar traits ($d$ is small) compete for the same resources, resulting in an energy cost for both.
- **Kernel:** $comp = \exp(-d^2 / 2\sigma^2_{comp})$.

### Mating
Agents with very similar traits can mate, resulting in an energy bonus (representing reproductive success or social cooperation).
- **Condition:** If $\exp(-d^2 / 2\sigma^2_{mate}) > 0.8$.

### Predation
If the trait difference is large enough, one agent acts as a predator and the other as prey.
- **Condition:** If $d > \text{threshold}$, the higher-trait agent gains energy while the lower-trait agent loses it.

## 5. System Configuration

The entire system is controlled by a `Config` object, which allows the user to define:
- Population size and world dimensions.
- The active movement strategy and its parameters (noise, friction, etc.).
- Ecological constants like interaction ranges and thresholds.

This configuration-driven approach makes the simulation highly extensible for future trait additions (e.g., life history speed, defense investment).
