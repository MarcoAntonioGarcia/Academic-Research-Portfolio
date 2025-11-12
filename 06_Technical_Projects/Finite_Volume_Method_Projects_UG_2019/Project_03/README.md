# Project 03: 1D Transient Diffusion (FVM)

## Academic Retrospective (November 2019)

This repository contains the work for an academic project developed in **November 2019** at the *Universidad de Guanajuato*. It is maintained as a **historical record** to exemplify the author's trajectory in research and computational problem-solving.

The project uses the **Finite Volume Method (FVM)** to solve the **one-dimensional thermal diffusion equation** in a **transient (time-dependent) state**.

---

### Technical Focus

The core objective was to implement and compare different time integration methods:

| Scheme | Stability | Temporal Parameter ($\Theta$) |
| :--- | :---: | :--- |
| **Implicit** | **Unconditionally Stable** | $\Theta = 1$ |
| **Explicit** | Conditionally Stable ($\Delta t$ restricted) | $\Theta = 0$ |

The analysis concluded that the **Implicit method is the most convenient** due to its inherent stability.

---