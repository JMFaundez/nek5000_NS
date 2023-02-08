# Momentum equation in Nek5000 

Nek routines to compute the different terms in the momentum equation in strong and weak formulation.

## NS

$$
\begin{align}
\frac{\partial u_i}{\partial t}   &= -\frac{\partial p}{\partial x_i} + \frac{1}{Re}\Delta u_i + f_i\ - u_j \frac{\partial u_i}{\partial x_j}\\
\frac{\partial u_i}{\partial x_i}&=0
\end{align}
$$

## Semi-discrete NS in Nek5000

$$
\begin{align}
\mathbf{M}\frac{\partial u_i}{\partial t}   &= \mathbf{D}_i^T p - \frac{1}{Re}\mathbf{K}u_i + \mathbf{M}f_i\ - \mathbf{C}u_i\\
\mathbf{D}_i u_i&=0
\end{align}
$$

- $\mathbf{C}$: Convection operator
- $\mathbf{M}$: Mass matrix
- $\mathbf{K}$: Stiffness matrix
- $\frac{\partial u_i}{\partial t}$ : Computed using BDF2, considering current solution `vx` and two previous time steps `vxlag`. Then scaled by the mass matrix
- $\mathbf{D}_i^T p$ : `opgradt(px,py,pz,pr)` (in `navier1.f`) where `pr` is defined in the pressure mesh and the outpouts `px`, `py` and `pz` are defined in the velocity mesh. 
- $\mathbf{K} u_i$ : `wlaplacian(lapu,u,diff,1)` (in `navier1.f`)
- $M f_i$ : `makeuf` will call user defined forces and put them in `BFX`, `BFY` and `BFZ` 
- The convective term is computed in nek as: `convop(convu,u_i)` and then scaled by the mass matrix (check `advab` routine in `navier1.f`)


## Time integration: Fractional step method

$$
\begin{align}
\mathbf{H}u_{i}^*&=\mathbf{D}_{i}^Tp^{n} + h_i^{n+1}\\
\frac{b_0}{\Delta t}\mathbf{D}_i\mathbf{M}^{-1}\mathbf{D}_{i}^T(p^{n+1}-p^{n})&= \mathbf{D}_i u_i^*\\
u_{i}^{n+1} &= u_i^* + \frac{\Delta t}{b_0}\mathbf{M}^{-1}\mathbf{D}_{i}^T(p^{n+1}-p^{n})
\end{align}
$$

- Helmoltz operator: $\mathbf{H}= \frac{b_0}{\Delta t}\mathbf{M} + \frac{1}{Re}\mathbf{K}$
- $h_i^{n+1}=-\sum_{j=1}^k\frac{b_j}{\Delta t} \mathbf{M}u_i^{n+1-j} - \sum_{j=1}^k a_j\mathbf{C}u_i^{n+1-j} + \mathbf{M}f_i^n$
