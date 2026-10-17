# Physics-Constrained Dynamic Identification and Noise Estimation for Underactuated Systems

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%20--%20R2026a-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![CasADi](https://img.shields.io/badge/CasADi-v3.5%20--%20v3.8+-brightgreen.svg)](https://web.casadi.org/)
[![IPOPT](https://img.shields.io/badge/Solver-IPOPT-orange.svg)](https://coin-or.github.io/Ipopt/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Official MATLAB/Simulink and CasADi implementation for the paper:  
**"Physics-Constrained Dynamic Identification and Noise Estimation for Underactuated Systems"**  
**Ahmed Asker** $^{1,2,*}$, **Omar Salah** $^{3}$, and **Haitham El-Hussieny** $^{1}$  
$^{1}$ *Department of Mechatronics and Robotics Engineering, Egypt-Japan University of Science and Technology (E-JUST), New Borg El-Arab City, Alexandria 21934, Egypt*  
$^{2}$ *Department of Production Engineering and Mechanical Design, Faculty of Engineering, Mansoura University, Mansoura 35516, Egypt*  
$^{3}$ *Mechatronics Engineering Department, Faculty of Engineering, Assiut University, Assiut 71516, Egypt*  
$^{*}$ *Corresponding author: ahmed.asker@ejust.edu.eg*

---

## 📌 Overview

Dynamic parameter identification of **underactuated mechanical systems** (systems with fewer independent control actuators than generalized degrees of freedom) is an inherently ill-posed problem. Because passive, unactuated joints receive **zero external actuation force or torque** ($\bm{\tau}_u \equiv \mathbf{0}$), their corresponding regression sub-problems are strictly homogeneous:

$$W_u(\bm{q}, \dot{\bm{q}}, \ddot{\bm{q}}) \bm{\Theta}_u \approx \mathbf{0}$$

Under measurement noise, encoder quantization, or numerical acceleration differentiation errors, unconstrained least-squares estimators trivially minimize algebraic residuals by driving passive parameters toward zero ($\bm{\Theta}_u \to \mathbf{0}$). This results in **numerical null-space collapse** ($>99\%$ loss of estimated mass and inertia), negative mass matrices violating positive definiteness ($D(\bm{q}) \nsucc 0$), and immediate divergence during forward dynamic integration.

While classical **Unconstrained Momentum-Filtered Least Squares (UMF-LS)** bypasses numerical acceleration differentiation through generalized momentum filtering, it remains an unconstrained linear transformation that does not prevent parameter collapse under noise.

This repository provides an end-to-end **Physics-Constrained Identification Framework** based on algorithmic differentiation over computational Directed Acyclic Graphs (DAGs) in **CasADi**:
1. **Automated Physical Constraint Extraction:** Automatic derivation of minimal linear symmetry ($A_{\mathrm{sym}} \bm{\Theta} = \mathbf{0}$) and passivity skew-symmetry ($A_{\mathrm{skew}} \bm{\Theta} = \mathbf{0}$) constraints via algorithmic differentiation and SVD/QR factorization, Sylvester principal minor positive-definiteness constraints ($\Delta_k(D(\bm{q})) \ge \epsilon_p > 0$), and sector-bounded physical parameter constraints ($\bm{\Theta} \in \mathcal{B}_{\mathrm{phys}}$).
2. **Joint Profile Maximum Likelihood Estimation (MLE):** Simultaneous identification of dynamic/friction parameters and statistical equation-error residual variance ($\hat{\sigma}^2$) to parameterize downstream Kalman filter process disturbance covariances ($Q_\tau = \hat{\sigma}^2 I$) and state observers.
3. **Normal-Equation Condensation:** Condensation of 40,000+ temporal observation rows into compact Gram matrix sufficient statistics ($H = W^T W$, $\bm{g} = W^T \bm{\Gamma}$), cutting interior-point solve times from $6.44\text{ s}$ to **$13\text{ ms}$** (**$488.7\times$ speedup**).
4. **Comprehensive Benchmarking on Simulation & Experimental Data:**
   - **3-DOF Double Inverted Pendulum on a Cart (DIPC):** High-fidelity multi-body Simscape simulation under controlled stochastic measurement noise ($\sigma_{\text{true}} = 0.05\text{ N}$).
   - **2-DOF Rotary Inverted Furuta Pendulum (RIFP):** Experimental hardware dataset ($20,000$ empirical data samples at $1\text{ kHz}$) evaluated under full-length ($20\text{ s}$) and truncated sliding observation horizons ($0.3\text{ s} - 2.0\text{ s}$).

---

## 🗂 Repository Structure

```
CasADi Underactuated System Identification/
├── README.md                           # Documentation, replication instructions, benchmark tables
├── LICENSE                             # MIT Open-Source License
├── .gitignore                          # MATLAB / LaTeX / Git build artifact filter
│
├── main_identification.m               # Master replication script (runs Benchmarks 1 & 2, generates tables & figures)
├── Generate_All_Figures.m              # Standalone publication figure generator (generates Figures 1-8 instantly)
├── run_dipc_benchmark_standalone.m     # Standalone runner for 3-DOF DIPC simulation benchmark
├── run_furuta_benchmark_standalone.m   # Standalone runner for 2-DOF RIFP experimental benchmark
│
├── core/                               # Algorithmic differentiation and kinematics engines
│   ├── casadi_dh_dyn.m                 # CasADi SX-based Denavit-Hartenberg dynamics, Jacobians, and passivity
│   ├── build_regressor.m               # Symbolic-to-numeric regressor generator & QR base parameter reduction
│   ├── extract_constraints.m           # Automated extraction of symmetry, skew-symmetry & positive-definiteness
│   ├── DH_dyn.m                        # General symbolic DH Lagrangian dynamics modeling
│   ├── Cart_sys.m                      # Symbolic 3-DOF DIPC modeling & Euler-Lagrange verification
│   └── Furuta_sys.m                    # Symbolic & CasADi 2-DOF RIFP modeling & verification
│
├── models/                             # Dynamic and contact/friction models
│   └── nonlinear_friction.m            # Regularized hyperbolic tangent symmetric & asymmetric friction models
│
├── identification/                     # Estimation and simulation engines
│   ├── opti_estimator.m                # High-performance condensed CasADi Opti estimator (MLE + physics constraints)
│   ├── benchmark_solvers.m             # Baseline solvers: Direct OLS, Decoupled OLS, Aguilar-Sierra Filtered LS
│   ├── run_dipc_benchmark.m            # 3-DOF DIPC identification pipeline, noise injection, and RK4 validation
│   └── simulate_furuta.m               # 2-DOF RIFP multi-step forward dynamic RK4 simulator
│
├── data/                               # Experimental and simulation datasets
│   ├── furutadataexp.mat               # 2-DOF RIFP experimental dataset (20,000 samples @ 1 kHz)
│   ├── Sim1.mat                        # 3-DOF DIPC simulation dataset under multi-frequency excitation
│   ├── Sim2.mat                        # Supplementary DIPC validation trajectory
│   ├── load_furuta_data.m              # Zero-phase FIR Nuttall filtering & central differencing
│   └── identification_results.mat      # Pre-computed benchmark workspace for instant figure reproduction
│
├── simulink/                           # Multi-body simulation models
│   ├── Model_under1.slx                # Underactuated 3-DOF DIPC Simscape Multibody model
│   ├── Model1.slx                      # Fully actuated reference model
│   ├── Model_under2.slx                # Alternative underactuated configuration model
│   └── run_dipc_simulation.m           # Script to configure parameters and run Simulink model to produce Sim1.mat
│
└── paper/                              # Paper manuscript, LaTeX sources, figures, and tables
    ├── main.tex                        # Full paper LaTeX source (IEEE Transactions style)
    ├── references.bib                  # BibTeX references
    ├── IEEEtran.cls                    # Document class file
    ├── main.pdf                        # Compiled publication PDF
    ├── figures/                        # All 8 publication figures (.pdf vector & .png 300 DPI)
    └── tables/                         # Auto-generated LaTeX benchmark comparison tables
```

---

## ⚙️ Prerequisites & Dependencies

- **MATLAB:** R2022b or later (tested up to R2026a).
- **Toolboxes:**
  - Optimization Toolbox (optional, if using MATLAB solvers)
  - Simulink & Simscape Multibody (only required if re-simulating `.slx` models from scratch)
- **CasADi:** Version 3.5.5 or later (tested on CasADi v3.8+).
  - Download from: [https://web.casadi.org/get/](https://web.casadi.org/get/)
  - Add CasADi to your MATLAB path:
    ```matlab
    addpath('C:\path\to\casadi-matlab');
    import casadi.*
    ```
- **Solvers:** IPOPT (bundled automatically with CasADi binary distributions).

---

## 🚀 Quick Start & How to Replicate All Results

### 1. Instant Reproduction of All Publication Figures (No Re-Optimization Required)
If you wish to immediately inspect and re-generate all 8 publication figures (vector `.pdf` and 300-DPI `.png`) using the pre-computed benchmark results stored in `data/identification_results.mat`:
```matlab
% In MATLAB command window:
Generate_All_Figures
```
This generates:
- `fig1_experimental_trajectories.pdf/.png`: Overview of RIFP 1 kHz trajectories ($q_1, q_2, \tau_1$).
- `fig2_torque_tracking_comparison.pdf/.png`: Actuated joint dynamic torque tracking and reconstruction residuals.
- `fig3_nullspace_collapse_demonstration.pdf/.png`: Normalized parameter recovery ratios and shrinking window stability ($T = 0.3\text{ s}$ to $2.0\text{ s}$).
- `fig4_forward_simulation_validation.pdf/.png`: Multi-step forward RK4 dynamic simulation validation driven by experimental torque $\tau_1(t)$.
- `fig5_physics_constraints_validation.pdf/.png`: Inertia matrix determinant $\det(D(\bm{q}(t))) > 0$ positive-definiteness certification.
- `fig6_friction_characteristic_curves.pdf/.png`: Identified symmetric and asymmetric regularized nonlinear joint friction curves.
- `fig_dipc_nullspace.pdf/.png`: Subsystem-grouped relative parameter estimation error across the 3-DOF DIPC under measurement noise ($\sigma = 0.05$).
- `fig_dipc_tracking.pdf/.png`: DIPC forward dynamic simulation trajectories ($d_1, \theta_2, \theta_3, F_1$).

---

### 2. Full End-to-End Pipeline Execution
To execute all identification algorithms, solvers, and simulations from scratch:
```matlab
% In MATLAB command window:
main_identification
```
This master script automatically:
1. Runs **Benchmark 1 (3-DOF DIPC)**:
   - Ingests `data/Sim1.mat` ($N = 2,000$ samples @ $1\text{ kHz}$).
   - Injects Gaussian measurement noise ($\sigma_{\text{true}} = 0.05\text{ N}$, $\sigma^2_{\text{true}} = 2.50 \times 10^{-3}\text{ N}^2$).
   - Solves Direct OLS (U-LS), Filtered LS (UMF-LS), and the proposed Physics-Constrained estimator (PC-LS).
   - Simulates forward dynamics (RK4) and exports `table4_dipc_comparison.tex`, `fig_dipc_nullspace`, and `fig_dipc_tracking`.
2. Runs **Benchmark 2 (2-DOF RIFP Experimental Hardware)**:
   - Ingests `data/furutadataexp.mat` ($N = 20,000$ samples @ $1\text{ kHz}$).
   - Computes zero-phase Nuttall FIR filtering and central difference derivatives.
   - Solves Decoupled LS (D-LS), Direct Unconstrained LS (U-LS), Direct PC-CasADi, Aguilar Filtered LS Baseline (UMF-LS), and Filter-Based Constrained LS (FBC-LS).
   - Solves asymmetric nonlinear friction identification.
   - Integrates forward RK4 dynamics over a $2.5\text{ s}$ validation window.
   - Saves all variables into `data/identification_results.mat`.
   - Exports `table1_parameters.tex`, `table2_performance_metrics.tex`, and Figures 1 through 6.

---

### 3. Running Individual Benchmarks Independently

#### A. 3-DOF Underactuated DIPC Benchmark Only
```matlab
run_dipc_benchmark_standalone
```

#### B. 2-DOF Furuta Pendulum Experimental Benchmark Only
```matlab
run_furuta_benchmark_standalone
```

#### C. Re-running the DIPC Multi-Body Simulation in Simulink
To re-generate the ground-truth simulation dataset (`data/Sim1.mat`) directly from the Simscape Multibody model under the multi-frequency excitation force $F_1(t) = 80 \sin(2\pi t) + 40 \cos(5\pi t)$:
```matlab
cd('simulink');
run_dipc_simulation;
cd('..');
```

#### D. Verifying Symbolic Mechanics and Denavit-Hartenberg Derivations
To symbolically verify the Lagrangian equations of motion, Denavit-Hartenberg kinematic chains, and Coriolis skew-symmetry:
```matlab
% 3-DOF DIPC Symbolic Verification:
cd('core');
Cart_sys;

% 2-DOF RIFP Symbolic & CasADi Verification:
Furuta_sys;
cd('..');
```

---

## 📊 Summary of Benchmark Results

### Benchmark 1: 3-DOF Underactuated DIPC under Measurement Noise ($\sigma_{\text{true}} = 0.05$, $\sigma^2_{\text{true}} = 2.50 \times 10^{-3}$)

| Subsystem | Parameter | Expression | Ground Truth | U-LS | UMF-LS (Aguilar) | Proposed PC-LS | Rel. Error |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Cart ($L_1$)** | $x_{1,1}$ | $m_1+m_2+m_3$ | $8.1000$ | $8.1031$ | $8.1045$ | **$8.0985$** | $0.02\%$ |
| | $x_{1,2}$ | $-(a_2 m_3 + lc_2 m_2)$ | $-1.0125$ | $-1.0125$ | $-1.0127$ | **$-1.0125$** | $0.00\%$ |
| | $x_{1,3}$ | $-lc_3 m_3$ | $-0.3375$ | $-0.3366$ | $-0.3363$ | **$-0.3375$** | $0.00\%$ |
| | $x_{1,4}$ | $-lc_3 m_3$ | $-0.3375$ | $-0.3372$ | $-0.3372$ | **$-0.3375$** | $0.00\%$ |
| | $x_{1,5}$ | $a_2 m_3 + lc_2 m_2$ | $1.0125$ | $1.0125$ | $1.0125$ | **$1.0125$** | $0.00\%$ |
| | $x_{1,6}$ | $lc_3 m_3$ | $0.3375$ | $0.3375$ | $0.3376$ | **$0.3375$** | $0.00\%$ |
| | $x_{1,7}$ | $2 lc_3 m_3$ | $0.6750$ | $0.6749$ | $0.6750$ | **$0.6750$** | $0.00\%$ |
| | $x_{1,8}$ | $lc_3 m_3$ | $0.3375$ | $0.3375$ | $0.3375$ | **$0.3375$** | $0.00\%$ |
| **Pendulum 1 ($L_2$)** | $x_{2,1}$ | $-(lc_2 m_2 + a_2 m_3)$ | $-1.0125$ | $0.0234^*$ | $0.0301^*$ | **$-1.0125$** | $0.00\%$ |
| | $x_{2,2}$ | $-lc_3 m_3$ | $-0.3375$ | $-0.0293^*$ | $-0.0382^*$ | **$-0.3375$** | $0.00\%$ |
| | $x_{2,3}$ | $I_{2zz}+I_{3zz}+a_2^2 m_3+lc_2^2 m_2+lc_3^2 m_3$ | $0.5631$ | $0.0004^*$ | $0.0000^*$ | **$0.5616$** | $0.26\%$ |
| | $x_{2,4}$ | $2 a_2 lc_3 m_3$ | $0.3375$ | $-0.0003^*$ | $-0.0004^*$ | **$0.3366$** | $0.26\%$ |
| | $x_{2,5}$ | $I_{3zz}+lc_3^2 m_3$ | $0.1128$ | $-0.0001^*$ | $-0.0002^*$ | **$0.1125$** | $0.25\%$ |
| | $x_{2,6}$ | $a_2 lc_3 m_3$ | $0.1688$ | $0.0001^*$ | $0.0000^*$ | **$0.1683$** | $0.26\%$ |
| | $x_{2,7}$ | $-a_2 lc_3 m_3$ | $-0.1688$ | $-0.0001^*$ | $0.0000^*$ | **$-0.1683$** | $0.26\%$ |
| | $x_{2,8}$ | $-2 a_2 lc_3 m_3$ | $-0.3375$ | $0.0001^*$ | $0.0002^*$ | **$-0.3366$** | $0.26\%$ |
| | $x_{2,9}$ | $g(lc_2 m_2 + a_2 m_3)$ | $-9.9292$ | $-0.0113^*$ | $-0.0044^*$ | **$-9.9055$** | $0.24\%$ |
| | $x_{2,10}$ | $g lc_3 m_3$ | $-3.3097$ | $0.0062^*$ | $0.0079^*$ | **$-3.3019$** | $0.24\%$ |
| **Pendulum 2 ($L_3$)** | $x_{3,1}$ | $-lc_3 m_3$ | $-0.3375$ | $-0.0065^*$ | $-0.0099^*$ | **$-0.3375$** | $0.00\%$ |
| | $x_{3,2}$ | $I_{3zz}+lc_3^2 m_3$ | $0.1128$ | $0.0001^*$ | $0.0000^*$ | **$0.1125$** | $0.25\%$ |
| | $x_{3,3}$ | $a_2 lc_3 m_3$ | $0.1688$ | $-0.0001^*$ | $-0.0001^*$ | **$0.1683$** | $0.26\%$ |
| | $x_{3,4}$ | $I_{3zz}+lc_3^2 m_3$ | $0.1128$ | $-0.0000^*$ | $-0.0000^*$ | **$0.1125$** | $0.25\%$ |
| | $x_{3,5}$ | $a_2 lc_3 m_3$ | $0.1688$ | $-0.0000^*$ | $-0.0001^*$ | **$0.1683$** | $0.26\%$ |
| | $x_{3,6}$ | $g lc_3 m_3$ | $-3.3097$ | $0.0009^*$ | $0.0010^*$ | **$-3.3019$** | $0.24\%$ |
| **Summary Metrics** | **Overall Param Error** | $\|\hat{\bm{\Theta}} - \bm{\Theta}^*\|_2 / \|\bm{\Theta}^*\|_2$ | --- | $80.01\%$ | $80.06\%$ | **$0.19\%$** | --- |
| | **Link 2 Norm** | $\|\bm{\Theta}_{L2}\|_2$ (True: $10.5498$) | $10.5498$ | $0.0397^*$ | $0.0494^*$ | **$10.5248$** | $0.24\%$ |
| | **Link 3 Norm** | $\|\bm{\Theta}_{L3}\|_2$ (True: $3.3393$) | $3.3393$ | $0.0066^*$ | $0.0099^*$ | **$3.3314$** | $0.24\%$ |
| | **Noise Variance $\hat{\sigma}^2$** | $\times 10^{-3}\text{ N}^2$ (True: $2.50$) | $2.5000$ | $2.5282$ | $2.5295$ | **$2.5339$** | $1.36\%$ |
| | **Simulation RMSE** | $d_1\text{ (m)} / \theta_2\text{ (rad)} / \theta_3\text{ (rad)}$ | --- | Diverged | Diverged | **$0.0195 / 0.1723 / 0.1801$** | --- |

*\*Indicates numerical null-space parameter collapse ($>99.5\%$ loss of parameter norm).*

> [!NOTE]
> **Unreduced Block-Diagonal Regressor Rationale:**  
> The formulation intentionally adopts an unreduced, subsystem-partitioned LIP parameterization $\bm{\Theta} \in \mathbb{R}^{24}$ and block-diagonal regressor $W = \text{blkdiag}(L_1, L_2, L_3) \in \mathbb{R}^{3N \times 24}$ rather than a manually reduced base-parameter set. Because $W$ is block-diagonal, unconstrained least-squares solvers (U-LS and UMF-LS) decouple across subsystems, failing to exploit shared physical parameters across links (such as $D_{12}=D_{21}$ and $D_{13}=D_{31}$). This setup rigorously demonstrates that automated algorithmic differentiation (AD) automatically discovers and enforces the necessary cross-subsystem coupling constraints $A_{\mathrm{sym}}\bm{\Theta} = \mathbf{0}$ and $A_{\mathrm{skew}}\bm{\Theta} = \mathbf{0}$ (isolated via SVD/QR factorization), completely resolving parameter identifiability without manual symbolic base-parameter reduction.

---

### Benchmark 2: 2-DOF Rotary Furuta Pendulum (RIFP Experimental Hardware Data)

| Parameter | Physical Description | Units | D-LS | U-LS | UMF-LS (Aguilar) | FBC-LS (Proposed) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| $\theta_1$ | Arm Base Inertia | $\text{kg}\cdot\text{m}^2$ | $0.0475$ | $0.0467$ | $0.0604$ | **$0.0604$** |
| $\theta_2$ | Pendulum Coupled Inertia | $\text{kg}\cdot\text{m}^2$ | $0.0144$ | $0.0135$ | $0.0148$ | **$0.0148$** |
| $\theta_3$ | Mutual Cross Inertia | $\text{kg}\cdot\text{m}^2$ | $0.0141$ | $0.0124$ | $0.0179$ | **$0.0179$** |
| $\theta_4$ | Pendulum Inertia | $\text{kg}\cdot\text{m}^2$ | $0.0000^*$ | $0.0075$ | $0.0126$ | **$0.0126$** |
| $\theta_5$ | Gravitational Moment ($m_2 g l_2$) | $\text{N}\cdot\text{m}$ | $0.0000^*$ | $0.3117$ | $0.4878$ | **$0.4878$** |
| $f_{v1}$ | Arm Viscous Friction | $\text{N}\cdot\text{m}\cdot\text{s/rad}$ | $0.0145$ | $0.0163$ | $0.0090$ | **$0.0090$** |
| $f_{v2}$ | Pendulum Viscous Friction | $\text{N}\cdot\text{m}\cdot\text{s/rad}$ | $0.0000^*$ | $0.0014$ | $0.0006$ | **$0.0006$** |
| $f_{c1}$ | Arm Coulomb Stiction | $\text{N}\cdot\text{m}$ | $0.0097$ | $0.0074$ | $0.0181$ | **$0.0181$** |
| $f_{c2}$ | Pendulum Coulomb Stiction | $\text{N}\cdot\text{m}$ | $0.0000^*$ | $0.0041$ | $0.0086$ | **$0.0086$** |

> [!NOTE]
> **Experimental Ground Truth:**  
> Link masses and CAD measurements provide kinematic completeness; derived base parameters $\theta_1, \dots, \theta_5$ and friction coefficients represent the calibrated experimental ground-truth benchmark identified in Aguilar-Sierra et al. (2014).

#### Directional Asymmetric Joint Friction Identification (RIFP)
When identified with directional friction asymmetry, the regularized friction formulation resolves into:
- **Actuated Arm ($q_1$):**
  - Positive rotation: $f_{c1}^+ = 0.0193\text{ N}\cdot\text{m}$, $f_{v1}^+ = 0.0082\text{ N}\cdot\text{m}\cdot\text{s/rad}$
  - Negative rotation: $f_{c1}^- = 0.0169\text{ N}\cdot\text{m}$, $f_{v1}^- = 0.0090\text{ N}\cdot\text{m}\cdot\text{s/rad}$
- **Passive Pendulum Bearing ($q_2$):**
  - Symmetric dissipation: $f_{c2} = 0.0086\text{ N}\cdot\text{m}$, $f_{v2} = 0.0006\text{ N}\cdot\text{m}\cdot\text{s/rad}$

This asymmetric formulation captures subtle directional variations while retaining smooth $C^\infty$ differentiability via $\tanh(\dot{q}_i/\epsilon_v)$ for interior-point solvers.

---

### Robustness Under Shrinking Observation Horizons ($T = 0.3\text{ s}$ to $2.0\text{ s}$)

| Window $N$ | Duration $T$ (s) | Method | $\min \det(D)\text{ (kg}^2\cdot\text{m}^4)$ | $\theta_4\text{ (kg}\cdot\text{m}^2)$ | $f_{v1}\text{ (N}\cdot\text{m}\cdot\text{s/rad)}$ | Sim RMSE $q_2$ (rad) |
| :---: | :---: | :--- | :---: | :---: | :---: | :---: |
| **$300$** | **$0.3$** | UMF-LS | $-1.17 \times 10^{-5}$ *(Indefinite)* | $-0.0028$ *(Negative)* | $-0.0312$ *(Energy Gen)* | **Diverged** |
| | | **FBC-LS (Proposed)** | $\mathbf{+1.00 \times 10^{-4}}$ | $\mathbf{+0.0039}$ | $\mathbf{+0.0000}$ | $\mathbf{0.0784}$ |
| **$500$** | **$0.5$** | UMF-LS | $-3.67 \times 10^{-5}$ *(Indefinite)* | $-0.0043$ *(Negative)* | $-0.0561$ *(Energy Gen)* | $0.9377\text{ }(53.7^\circ)$ |
| | | **FBC-LS (Proposed)** | $\mathbf{+1.00 \times 10^{-4}}$ | $\mathbf{+0.0044}$ | $\mathbf{+0.0000}$ | $\mathbf{0.0617}\text{ }(3.5^\circ)$ |
| **$800$** | **$0.8$** | UMF-LS | $+1.68 \times 10^{-4}$ | $+0.0089$ | $-0.0124$ *(Energy Gen)* | $0.2140$ |
| | | **FBC-LS (Proposed)** | $\mathbf{+2.93 \times 10^{-4}}$ | $\mathbf{+0.0112}$ | $\mathbf{+0.0042}$ | $\mathbf{0.0381}$ |
| **$1000$** | **$1.0$** | UMF-LS | $+1.85 \times 10^{-4}$ | $+0.0094$ | $-0.0068$ *(Energy Gen)* | $0.1852$ |
| | | **FBC-LS (Proposed)** | $\mathbf{+3.45 \times 10^{-4}}$ | $\mathbf{+0.0118}$ | $\mathbf{+0.0065}$ | $\mathbf{0.0246}$ |
| **$2000$** | **$2.0$** | UMF-LS | $+4.06 \times 10^{-4}$ | $+0.0121$ | $+0.0082$ | $0.0152$ |
| | | **FBC-LS (Proposed)** | $\mathbf{+4.06 \times 10^{-4}}$ | $\mathbf{+0.0124}$ | $\mathbf{+0.0088}$ | $\mathbf{0.0114}$ |

---

## 🔬 Mathematical Formulations

### 1. Partitioned Underactuated Dynamics
Mechanical systems with $n$ generalized coordinates and $m < n$ independent control actuators partition into actuated ($\bm{q}_a \in \mathbb{R}^m$) and unactuated ($\bm{q}_u \in \mathbb{R}^{n-m}$) coordinates:

$$\begin{bmatrix} D_{aa}(\bm{q}) & D_{au}(\bm{q}) \\ D_{ua}(\bm{q}) & D_{uu}(\bm{q}) \end{bmatrix} \begin{bmatrix} \ddot{\bm{q}}_a \\ \ddot{\bm{q}}_u \end{bmatrix} + \begin{bmatrix} \bm{C}_{aa}\dot{\bm{q}}_a + \bm{C}_{au}\dot{\bm{q}}_u \\ \bm{C}_{ua}\dot{\bm{q}}_a + \bm{C}_{uu}\dot{\bm{q}}_u \end{bmatrix} + \begin{bmatrix} \bm{g}_a(\bm{q}) \\ \bm{g}_u(\bm{q}) \end{bmatrix} + \begin{bmatrix} \bm{\tau}_{fa}(\dot{\bm{q}}_a) \\ \bm{\tau}_{fu}(\dot{\bm{q}}_u) \end{bmatrix} = \begin{bmatrix} \bm{\tau}_a \\ \mathbf{0} \end{bmatrix}$$

Because unactuated coordinates receive zero generalized control forces ($\bm{\tau}_u \equiv \mathbf{0}$), their identification equation is strictly homogeneous ($W_u \bm{\Theta} \approx \mathbf{0}$).

### 2. Algorithmic Differentiation & Minimal Constraint Extraction via SVD/QR
To eliminate polynomial explosion in symbolic computer algebra systems, the framework builds computational Directed Acyclic Graphs (DAGs) in CasADi and extracts exact constraint matrices automatically:
- **Linear Inertia Symmetry ($A_{\mathrm{sym}} \bm{\Theta} = \mathbf{0}$):** Evaluates row gradients:
  $$\bm{a}_{\mathrm{sym}, ij}^T = \frac{\partial \left( D_{ij}(\bm{q}) - D_{ji}(\bm{q}) \right)}{\partial \bm{\Theta}}$$
  across independent configurations on the state manifold and extracts the minimal linearly independent row space via Singular Value Decomposition (SVD) or rank-revealing QR factorization.
- **Passivity Skew-Symmetry ($A_{\mathrm{skew}} \bm{\Theta} = \mathbf{0}$):** The fundamental physical passivity condition $\dot{\bm{q}}^T (\dot{D} - 2C) \dot{\bm{q}} = 0$ requires skew-symmetry $N(\bm{q}, \dot{\bm{q}}) = \dot{D} - 2C = -N^T$. The minimal constraint matrix $A_{\mathrm{skew}}$ is extracted by evaluating the row gradients of skew-symmetric pairs across test configurations and isolating the linearly independent basis via SVD/QR.
- **Sylvester Principal-Minor Positive Definiteness ($\Delta_k(D(\bm{q})) \ge \epsilon_p$):** Strict positive definiteness $D(\bm{q}) \succ 0$ is enforced configuration-wise across the reachable operational envelope $\mathcal{Q}_{\text{reach}}$, ensuring well-conditioned numerical inversion throughout forward simulation and MPC without non-convex 10-parameter pseudo-inertia reconstructions.

### 3. Profile Maximum Likelihood Estimation & Equation-Error Residuals
Assuming additive Gaussian measurement noise on the canonical equation-error observation model:

$$\bm{\Gamma}_k = W_k \bm{\Theta} + \bm{\varepsilon}_k, \quad \bm{\varepsilon}_k \sim \mathcal{N}(\mathbf{0}, \sigma^2 I_{d \times d})$$

The negative log-likelihood across all $N$ temporal observations ($M = d \cdot N$ total scalar measurements) is:

$$\mathcal{L}(\bm{\Theta}, \sigma^2) = \frac{d N}{2} \ln(2\pi) + \frac{d N}{2} \ln(\sigma^2) + \frac{1}{2\sigma^2} \| W\bm{\Theta} - \bm{\Gamma} \|_2^2$$

Concentrating out $\sigma^2$ analytically via its optimal variance estimate $\hat{\sigma}^2(\bm{\Theta}) = \frac{1}{d N} \| W\bm{\Theta} - \bm{\Gamma} \|_2^2$ produces the profile log-likelihood objective:

$$\min_{\bm{\Theta}} \quad J(\bm{\Theta}) = \frac{d N}{2} \ln \left( \| W\bm{\Theta} - \bm{\Gamma} \|_2^2 \right) \equiv \frac{d N}{2} \ln \left( \bm{\Theta}^T H \bm{\Theta} - 2 \bm{g}^T \bm{\Theta} + c \right)$$

where $H = W^T W$, $\bm{g} = W^T \bm{\Gamma}$, and $c = \bm{\Gamma}^T \bm{\Gamma}$ are pre-computed BLAS-3 sufficient statistics.

### 4. Complete Optimization Problem with Physical Admissible Set $\mathcal{B}_{\mathrm{phys}}$
$$\begin{aligned}
\min_{\bm{\Theta}} \quad & \frac{d N}{2} \ln \left( \bm{\Theta}^T H \bm{\Theta} - 2 \bm{g}^T \bm{\Theta} + c \right) \\
\text{s.t.} \quad & A_{\mathrm{sym}} \bm{\Theta} = \mathbf{0} \\
& A_{\mathrm{skew}} \bm{\Theta} = \mathbf{0} \\
& \Delta_k(D(\bm{q}^{(j)}; \bm{\Theta})) \ge \epsilon_p, \quad k = 1, \dots, n, \quad j = 1, \dots, N_c \\
& \bm{\Theta} \in \mathcal{B}_{\mathrm{phys}}
\end{aligned}$$

where the physically admissible parameter set is defined as:

$$\mathcal{B}_{\mathrm{phys}} = \left\{ \bm{\Theta} : S_{\text{pos}}\bm{\Theta}_{\text{dyn}} \ge \mathbf{0}, \quad \bm{\Theta}_{\text{fric}} \ge \mathbf{0}, \quad \bm{\Theta}_{\text{sector}} \le \mathbf{0} \right\}$$

Here, $S_{\text{pos}}$ selects strictly positive physical parameters (such as cart mass, link inertias, and friction coefficients), while admitting negative sector bounds for gravitational and geometric cross-terms (such as $x_{2,9}, x_{2,10}, x_{3,6} \le 0$ on the DIPC).

### 5. Downstream Dynamic Observer Parameterization & Covariance Calibration
The estimated residual variance $\hat{\sigma}^2$ directly parameterizes state estimators and Kalman filters:
- **Dynamic Process Disturbance Covariance:** In robotics state observers (EKF, UKF, MHE), generalized torque disturbances are parameterized by $Q_\tau = \hat{\sigma}^2 I_{d \times d}$ (or torque measurement noise covariance $R_\tau = \hat{\sigma}^2 I_{d \times d}$ when direct force/torque sensing is present).
- **Parameter Estimation Covariance:** Calibrated parameter uncertainty ellipsoids are obtained directly from the inverted Hessian at the constrained optimum:
  $$\text{Cov}(\hat{\bm{\Theta}}) \approx \left( \nabla^2 J(\hat{\bm{\Theta}}) \right)^{-1} \approx \hat{\sigma}^2 H^{-1}$$
  or projected onto the null space of active linear equality constraints $\mathcal{N}([A_{\mathrm{sym}}^T, A_{\mathrm{skew}}^T]^T)$.

---

## 📖 Citation

If you use this codebase or benchmark data in your research, please cite our paper and repository:

```bibtex
@article{asker2026physics,
  title={Physics-Constrained Dynamic Identification and Noise Estimation for Underactuated Systems},
  author={Asker, Ahmed and Salah, Omar and El-Hussieny, Haitham},
  journal={},
  year={2026},
  note={Under review}
}

@misc{AskerRepo2026,
  author       = {Asker, Ahmed},
  title        = {{CasADi Underactuated System Identification: Physics-Constrained Dynamic Identification and Noise Estimation for Underactuated Systems}},
  year         = {2026},
  howpublished = {\url{https://github.com/Ahmed-Askar-GB/CasADi-Underactuated-System-Identification}},
  note         = {{GitHub} repository}
}
```

---

## 📄 License
This repository is released under the [MIT License](LICENSE).
