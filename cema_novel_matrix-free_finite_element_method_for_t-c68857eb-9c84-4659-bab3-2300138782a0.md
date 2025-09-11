# A Novel Matrix-Free Finite Element Method for Time-Harmonic Maxwell's Equations 

Yi-Yao Wang, Graduate Student Member, IEEE, Qiwei Zhan ${ }^{\text {® }}$, Member, IEEE, Haoqiang Feng, Yin-Da Wang ${ }^{\text {® }}$, Student Member, IEEE, Bozhao Sun, and Wen-Yan Yin ${ }^{\text {® }}$, Fellow, IEEE


#### Abstract

A frequency-domain finite element method (FEM) with the low-storage matrix-free feature is proposed for efficient analysis of electromagnetic fields. As opposed to conventional frequency-domain FEM which requires the matrix assembly before solving, the proposed matrix-free algorithm avoids the assembly and storage of the global matrix. In this approach, the global sparse matrix-vector (SpMV) multiplication is decomposed into element-wise matrix-vector (MV) multiplications. Additionally, the sum factorization technique is applied on tensorial basis functions to reduce the complexities of local MV multiplications. The numerical results demonstrate the superiority of our algorithm over traditional FEM in terms of both memory and time consumption. Besides, the improvement is more profound when higher-order basis functions are considered. A speedup of more than 10 against matrix-assembled solvers is observed. Given that memory transfer is bounded more than computation resources in modern supercomputer architectures, our algorithm is more friendly to high-performance computing platforms.


Index Terms-Geometric multigrid, matrix-free, sum factorization, tensorial basis functions, vectorization.

## I. Introduction

HIGH-PERFORMANCE computational electromagnetic methods (CEMs) are crucial for advanced applications such as antenna design, microwave device optimization, electromagnetic protection capability assessment, etc. With the ascending frequency of interest, the progressively challenging problems necessitate advanced CEM to exhibit the following features: scalability, adaptivity, efficiency, and stability [1]. Widely adopted numerical methods include the method of moments [2], finite-difference time-domain [3], finite element method (FEM) [4], and discontinuous Galerkin FEM (DG-FEM) [5], [6]. Of the various methods above, FEM is particularly advantageous due to its capability of modeling intricate geometries, adaptability to solve both electrical small

[^0]and large problems [7], and convenient coupling with multiphysical fields [8]. Consequently, FEM has been adopted as the method of choice in this study.

In [5] and [6], we have proposed an efficient long-time stable DG time-domain (DGTD) method with explicit timestepping. By combining the Gauss-Lobatto basis function and quadrature points on tensor-product elements [9], we conducted a fast and low-storage time-domain algorithm, avoiding the expensive solution of linear systems. In the case of frequency-domain methods, however, the assembly and solution of a large linear system are usually inevitable. Emerging as a popular topic in the last decade, the matrix-free approach enables the possibility of avoiding assembling a global matrix even for nontemporal problems [10]. The matrix-free approach exploits the nature that the global matrix only appears in products with other vectors if iterative solvers are applied. Memory-bounded [11] global SpMV multiplication is decomposed as element-wise MV computations. This technique has been widely adopted by some projects [12], [13], [14], [15] and applied to thermal analysis [16], fluid-dynamics [13], [14], [17], structural mechanics [18], etc. Kronbichler and Kormann [19] proposed a matrix-free DG method that reaches up to $60 \%$ of the CPU arithmetic peak. They also extend matrix-free methods to GPU architectures [13], implicit Runge-Kutta methods [20], and multigrid preconditioner on locally-refined mesh [21].

Though requiring lower memory, matrix-free methods must compute the element-wise product of matrices and vectors in each iteration. The sum factorization technique, first proposed in the spectral element method, reduces the complexity of local matrix assembly and MV multiplications [9], [22], [23]. Sum factorization is typically implemented on hexahedral and quadrilateral elements with tensor-product basis functions, while its application on a simplicial mesh is still an open area. Moxey et al. [24], [25] demonstrate the sum factorization on tetrahedron, pyramid, and prism by mapping their basis functions to the collapsed hexahedron in Nektar ++ and implement matrix-free method on a simplicial mesh. Sun et al. [26] compare the performance of matrix-free methods between hexahedral/quadrilateral and tetrahedral/triangular mesh on positive Helmholtz and elasticity operators. However, throughput still has a significant drop on simplicial mesh in both [25] and [26]. Overall, matrix-free methods have shown promising performance over traditional FEM for high-order methods, because it's more arithmetic-intensive and requires lower memory.

Despite its potential, the matrix-free method has seldom been employed for electromagnetic problems. In this article, a matrix-free FEM suitable for high-order methods is proposed for the frequency-domain Maxwell's wave equation. It's necessary to distinguish our work from other algorithms also known as matrix-free in the CEM community, such as time-domain FEM [5], time-domain based frequency-domain method [27], or multilevel fast multipole algorithm (MLFMA) [28]. Since the matrix-free method on simplicial mesh is less efficient, this work will rely on quadrilateral and hexahedral discretization in 2-D and 3-D, respectively.

This article is outlined as follows. In Section II, we review the discretized weak form of frequency-domain FEM and analyze the complexity of matrix assembly with respect to basis order. Section III proposes the element-based matrixfree method, the evaluation of local MV product on reference cell, and the sum factorization technique. In Section IV, a geometric multigrid preconditioner available under the matrix-free framework is suggested. Then we demonstrate the superior runtime and memory performance of matrix-free frequencydomain FEM on different numerical examples in Section V. Finally, discussions and probable improvements are given in Section VI.

## II. Traditional FEM in Frequency-Domain

FEM is effective for the near-field solution of electromagnetic problems with complex geometries. Some emerging FEM topics within the field of electromagnetics include hybridizable discontinuous Galerkin method [29], [30] and discontinuous Petrov-Galerkin method [31]. In this section, we review the process of conventional matrix-assembled FEM and analyze the complexity of matrix assembly.

## A. FEM Formulation of Maxwell's Equations

1) Strong Form: Frequency-domain solvers are preferred for devices with fine multiscale structures where time-domain methods may require unduly short time steps. Suppose harmonically oscillating functions and eliminate the magnetic terms, we have the time-harmonic Maxwell's wave equation

$$
\begin{equation*}
\nabla \times\left(\frac{1}{\mu_{r}} \nabla \times \boldsymbol{E}\right)-k_{0}^{2} \epsilon_{r} \boldsymbol{E}=-j k_{0} Z_{0} \boldsymbol{J}, \text { in } \Omega . \tag{1}
\end{equation*}
$$

Common boundary conditions such as perfect electric conductor (PEC) or absorption boundary condition (ABC) are applied on boundaries $\Gamma_{1}$ and $\Gamma_{2}$ respectively,

$$
\begin{align*}
\boldsymbol{n} \times \boldsymbol{E} & =0, \text { on } \Gamma_{1}  \tag{2}\\
\boldsymbol{n} \times\left(\frac{1}{\mu_{r}} \nabla \times \boldsymbol{E}^{s}\right)+j k_{0} \boldsymbol{n} \times\left(\boldsymbol{n} \times \boldsymbol{E}^{s}\right) & =0, \text { on } \Gamma_{2} \tag{3}
\end{align*}
$$

where $\boldsymbol{E}$ is the total electric field, $\boldsymbol{E}^{s}$ is the scattered electric field going outward, $\boldsymbol{J}$ represents current source, $k_{0}$ is wavenumber in vacuum, $\mu_{r}$ and $\epsilon_{r}$ are relative permeability and permittivity, respectively. $Z_{0}$ is the vacuum impedance, and $\boldsymbol{n}$ is the normal vector on boundaries.

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-02.jpg?height=264&width=857&top_left_y=184&top_left_x=1092)
Fig. 1. Vector nodal basis functions on a quadrilateral element.

2) Weak Form: The exact solution of (1)-(3) is difficult to find analytically. Suppose $\boldsymbol{E} \in \boldsymbol{U}$ where $\boldsymbol{U}$ is the trial space. (1) is tested by another function $\boldsymbol{W}$ in the test space $\boldsymbol{V}$. Applying integral by parts and Green's theorem, we have

$$
\begin{align*}
& \left(\nabla \times \boldsymbol{W}, \frac{1}{\mu_{r}} \nabla \times \boldsymbol{E}\right)_{\Omega}-k_{0}^{2} \epsilon_{r}(\boldsymbol{W}, \boldsymbol{E})_{\Omega} \\
& \quad-\left\langle\boldsymbol{W}, j \frac{k_{0}}{\eta_{r}} \boldsymbol{n} \times \boldsymbol{n} \times \boldsymbol{E}\right\rangle_{\Gamma_{2}}=-\left(\boldsymbol{W}, j k_{0} Z_{0} \boldsymbol{J}\right)_{\Omega} \tag{4}
\end{align*}
$$

where $(\cdot, \cdot)_{\Omega}$ is volume integral in domain $\Omega$ and $\left\langle\cdot, \cdot>_{\Gamma}\right.$ is surface integral on $\Gamma$. $\eta_{r}=\left(\mu_{r} / \epsilon_{r}\right)^{1 / 2}$ is the relative impedance.
3) Discretized Weak Form: The FEM finds the solution of $\boldsymbol{E}$ on $\Omega_{K}$, the discretized domain of $\Omega$. A field within an element is interpolated by a set of polynomial basis functions. In the Galerkin method, the function space $\boldsymbol{V}$ is chosen the same as $\boldsymbol{U}$. The discretized electric field at $\boldsymbol{r}$ is

$$
\begin{align*}
\boldsymbol{E}_{h}(\boldsymbol{r}) & =\sum_{i=1}^{|\Phi|} E_{i} \boldsymbol{\Phi}_{i}(\boldsymbol{r}) \\
& =\sum_{i=1}^{|\Phi|} E_{i} \Phi_{i}(\boldsymbol{r}) \hat{\boldsymbol{e}}_{i} \tag{5}
\end{align*}
$$

where $\Phi_{i}(\boldsymbol{x})$ and $\hat{\boldsymbol{e}}_{i}$ are basis functions and the unit Cartesian direction vector corresponding to the $i$ th freedom, respectively. $|\Phi|$ defines the number of global DoF. Fig. 1 depicts the complete nodal vector bases on a quadrilateral element.

## B. Linear System Assembly of discretized Weak Form

An important concept in FEM is assembly. It involves the aggregation of local stiffness matrices and right-hand side vectors into a global linear system $\boldsymbol{A x}=\boldsymbol{b}$, where

$$
\begin{align*}
\boldsymbol{A}_{i j}= & \int_{\Omega_{K}}\left[\nabla \times \boldsymbol{\Phi}_{i} \cdot \frac{1}{\mu_{r}} \nabla \times \boldsymbol{\Phi}_{j}-k_{0}^{2} \epsilon_{r} \boldsymbol{\Phi}_{i} \cdot \boldsymbol{\Phi}_{j}\right] d x \\
& +\int_{\partial \Omega_{K}}\left[j \frac{k_{0}}{\eta_{r}}\left(\boldsymbol{n} \times \boldsymbol{\Phi}_{i}\right) \cdot\left(\boldsymbol{n} \times \boldsymbol{\Phi}_{j}\right)\right] d s \tag{6}
\end{align*}
$$

$\boldsymbol{x}$ is the unknowns and $\boldsymbol{b}$ stands for the contribution of excitation. In high-order methods, the assembly of $\boldsymbol{A}$ is a bottleneck because both assembly time and storage are substantially higher than in linear methods. For the convenience of later deductions in Section III, we will consider a Helmholtz equation which is the scalar form of (1)

$$
\begin{equation*}
-\nabla^{2} u-k^{2} u=f \quad \text { in } \Omega \tag{7}
\end{equation*}
$$

Different to [25] and [26] where they performed tests on matrix-free positive-definite Helmholtz operators, (7) is indefinite. $\boldsymbol{A}$ is stored in compressed sparse row (CSR) format [11].

TABLE I
Complexities of Local Matrix Assembly, mV Evaluation, and Matrix Storage in Traditional FEM
|  | conventional | optimized |
| :---: | :---: | :---: |
| Assembly | $O\left(P^{3 d}\right)$ | $O\left(P^{2 d+1}\right)$ |
| Matrix-vector | $O\left(P^{2 d}\right)$ | $O\left(P^{d+1}\right)$ |
| Storage | $O\left(P^{2 d}\right)$ | $O\left(P^{d}\right)$ |


A restriction operator $\boldsymbol{R}_{K}$ is defined that maps the global DoF into local ones in element $K . \boldsymbol{R}_{K}$ is a rectangular matrix with either 0 or 1 entries and $\boldsymbol{R}_{K}^{T}$ projects the local DoF back to global DoF. Let $\mathcal{B}$ and $\mathcal{Q}$ be the number of basis functions and quadrature points in an element, respectively. Decompose the global matrix to its element-wise counterparts $\boldsymbol{A}_{K} \in \mathbb{R}^{\mathcal{B} \times \mathcal{B}}$, we have

$$
\begin{align*}
\boldsymbol{A} & =\sum_{K=1}^{|K|} \boldsymbol{R}_{K}^{T} \boldsymbol{A}_{K} \boldsymbol{R}_{K} \\
\boldsymbol{A}_{K} & =\boldsymbol{W}_{K}-k^{2} \boldsymbol{M}_{K} \tag{8}
\end{align*}
$$

and

$$
\begin{align*}
\left(\boldsymbol{W}_{K}\right)_{i j} & =\int_{\Omega_{K}} \nabla \Phi_{i}(\boldsymbol{x}) \cdot \nabla \Phi_{j}(\boldsymbol{x}) d x \\
& =\sum_{q=1}^{\mathcal{Q}} \nabla \Phi_{i}\left(\boldsymbol{x}_{q}\right) \cdot \nabla \Phi_{j}\left(\boldsymbol{x}_{q}\right) w_{q} \\
\left(\boldsymbol{M}_{K}\right)_{i j} & =\int_{\Omega_{K}} \Phi_{i}(\boldsymbol{x}) \Phi_{j}(\boldsymbol{x}) d x \\
& =\sum_{q=1}^{\mathcal{Q}} \Phi_{i}\left(\boldsymbol{x}_{q}\right) \Phi_{j}\left(\boldsymbol{x}_{q}\right) w_{q} \tag{9}
\end{align*}
$$

where $|K|$ represents the number of cells, $\boldsymbol{W}_{K}$ and $\boldsymbol{M}_{K}$ are the derivative and mass matrices, respectively, and $w_{q}$ is the quadrature weight at $q$ th quadrature point.

It is obvious from (9) that the assembly of $\boldsymbol{A}_{K}$ goes through three loops: one for quadrature points, one for test functions, and one for trial functions. For scalar nodal basis functions with exact integration [5], [32]

$$
\begin{equation*}
\mathcal{B}=\mathcal{Q}=P^{d}, \quad d=1,2,3 \tag{10}
\end{equation*}
$$

where $P$ is the number of basis functions in one dimension and $d$ represents the dimension. For Lagrange basis functions with $p$ th order polynomials, $P=p+1$. The computational complexity of local matrix assembly is $O\left(P^{3 d}\right)$ and the storage complexities of a local matrix is $O\left(P^{2 d}\right)$. Consider an example where $p$ rises from 1 to 2 , the calculation cost of one local matrix will exhibit a growth of over 11 -fold in 2-D and 38fold in 3-D. The rapid increase in both time and storage inevitably undermines the advantage of high-order methods [6]. The situation is even more severe for Nédélec basis functions because $\mathcal{B}=3 p(p+1)^{2}$ in hexahedral elements [33]. If matrix assembly is optimized with sum factorization, as will be shown in Section III, then the complexity could be reduced to $O\left(P^{2 d+1}\right)$. Computational complexities for different element-wise operations and storage complexity concerning $P$ are provided in Table I.

Equivalently, we can rewrite the assembly of $\boldsymbol{W}_{K}$ and $\boldsymbol{M}_{K}$ as products of matrices [34], [35]

$$
\begin{align*}
\boldsymbol{W}_{K} & =\boldsymbol{G}^{T} \boldsymbol{D}_{d} \boldsymbol{G} \\
\boldsymbol{M}_{K} & =\boldsymbol{B}^{T} \boldsymbol{D} \boldsymbol{B} \tag{11}
\end{align*}
$$

where

$$
\begin{align*}
\boldsymbol{D} & =\operatorname{diag}\left(w_{1}, w_{2}, \ldots, w_{\mathcal{Q}}\right) \\
\boldsymbol{D}_{d} & =\boldsymbol{D} \otimes \boldsymbol{I}_{d} \tag{12}
\end{align*}
$$

are diagonal matrices of 1 or $d$ copies of quadrature weight times Jacobian determinant, and

$$
\begin{align*}
\boldsymbol{B}_{q i} & =\Phi_{i}\left(\boldsymbol{x}_{q}\right) \\
\boldsymbol{G} & =\left[\boldsymbol{G}_{1}^{T} \cdots \boldsymbol{G}_{\mathcal{Q}}^{T}\right]^{T} \\
\left(\boldsymbol{G}_{q}\right)_{j i} & =\left.\frac{\partial \Phi_{i}}{\partial x_{j}}\right|_{\boldsymbol{x}_{q}}, \quad q \in[1, \mathcal{Q}], \quad i \in[1, \mathcal{B}], \quad j \in[1, d] \tag{13}
\end{align*}
$$

are matrices representing basis function values or gradient values at the $q$ th quadrature point, respectively. The procedure of a prototypical SpMV computation is shown in Algorithm 1.

## C. Solution of the Linear System

Linear system $\boldsymbol{A x}=\boldsymbol{b}$ can be solved with either direct methods or iterative methods. Sparse direct solvers invert $\boldsymbol{A}$ with a LU decomposition. They quickly become infeasible for large-scale 3-D problems when the number of DoF exceeds several million. Iterative solvers, on the other hand, seek an approximate solution. Instead of inverting $\boldsymbol{A}$ directly, only MV multiplications are required. For example, Krylov subspace methods find a solution $\boldsymbol{x}=q_{m-1}(\boldsymbol{A}) \boldsymbol{b} \approx \boldsymbol{A}^{-1} \boldsymbol{b}$, in which $q_{m-1}$ is a polynomial of degree $m-1$ [36]. This makes iterative methods more suitable for large-scale problems and parallelization on high-performance computing platforms. We will demonstrate how this key property is leveraged to avoid constructing and storing the global matrix in Section III.

```
Algorithm 1 Matrix-Assembled SpMV Multiplication
    input : $\boldsymbol{v} \in \mathbb{R}^{n}$
    output : $\boldsymbol{u} \in \mathbb{R}^{n}$
    initialize: Sparsity pattern of $\boldsymbol{A}$
    for $K=1$ to $|K|$ do
        initialize: $\boldsymbol{R}_{K}, \boldsymbol{B}, \boldsymbol{G}, \boldsymbol{D}$, and $\boldsymbol{D}_{d}$
        $\boldsymbol{A}_{K} \leftarrow \boldsymbol{G}^{T} \boldsymbol{D}_{d} \boldsymbol{G}-k^{2} \boldsymbol{B}^{T} \boldsymbol{D} \boldsymbol{B}$;
        $\boldsymbol{A} \leftarrow \boldsymbol{A}+\boldsymbol{R}_{K}^{T} \boldsymbol{A}_{K} \boldsymbol{R}_{K}$
    $\boldsymbol{u} \leftarrow \boldsymbol{A v}$
```


## III. Matrix-Free Method for Frequency-Domain Electromagnetic FEM

High-order finite elements entail a higher computational and storage complexity to construct and store $\boldsymbol{A}$. Besides, SpMV is limited by the memory bandwidth which grows lower than arithmetic capabilities and the number of computing cores. A SpMV has been shown to achieve less than $10 \%$ of peak on 2010's hardware [11], and with the increase in

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-04.jpg?height=644&width=803&top_left_y=178&top_left_x=203)
Fig. 2. Process of an SpMV multiplication decomposition. (a) Mesh consists of 2 quadrilateral elements with first-order basis functions. (b) Element-wise decomposition of the global SpMV multiplication. The restriction and prolongation operators are omitted.

arithmetic performance during the last decade the achievable peak is around $0.5 \%-2 \%$. These shortcomings limit the wide application of cutting-edge high-order methods [6].

## A. Element-Wise Decomposition of SpMV

Instead of factorized directly, the global matrix only appears as the product with other vectors in iterative solvers. For example, the solution of Richardson iteration is adapted as [36]

$$
\begin{align*}
\boldsymbol{r}_{k} & :=\boldsymbol{b}-\boldsymbol{A} \boldsymbol{x}_{k} \\
\boldsymbol{x}_{k+1} & :=\boldsymbol{x}_{k}+\boldsymbol{A} \boldsymbol{r}_{k} \tag{14}
\end{align*}
$$

where $\boldsymbol{r}_{k}$ is residual at previous step. In traditional FEM, SpMV products require the constructed global matrix. However, as will be shown in the numerical results, storage of the global matrix may take up more than half of a program's memory. Alternatively, by decomposing both the global matrix and input vector as their element-wise counterparts, the SpMV multiplication is computed by [10]

$$
\begin{align*}
\boldsymbol{A} \boldsymbol{v} & =\left(\sum_{K=1}^{|K|} \boldsymbol{R}^{T} \boldsymbol{A}_{K} \boldsymbol{R}_{K}\right) \boldsymbol{v} \\
& =\sum_{K=1}^{|K|} \boldsymbol{R}_{K}^{T} \boldsymbol{A}_{K} \boldsymbol{R}_{K} \boldsymbol{v} \\
& =\sum_{K=1}^{|K|} \boldsymbol{R}_{K}^{T} \boldsymbol{A}_{K} \boldsymbol{v}_{K} . \tag{15}
\end{align*}
$$

Fig. 2 illustrates this process on a quadrilateral mesh with 6 basis functions.

With (15), multiplication of $\boldsymbol{A} \boldsymbol{v}$ does not need explicit storage of $\boldsymbol{A}$ if local matrices $\boldsymbol{A}_{K}$ are computed on the fly. Two considerations must be taken here. Firstly, the matrix-matrix product in traditional is avoided by evaluating $\boldsymbol{A}_{K} \boldsymbol{v}_{K}$ from right to left, as shown in Algorithm 2. Compared to Algorithm 1, the complexity is reduced from $O\left(P^{3 d}\right)$

```
Algorithm 2 Matrix-Free MV Multiplication
    input : $\boldsymbol{v} \in \mathbb{R}^{n}$
    output: $\boldsymbol{u} \in \mathbb{R}^{n}$
    for $K=1$ to $|K|$ do
        initialize: $\boldsymbol{R}_{K}, \boldsymbol{B}, \boldsymbol{G}, \boldsymbol{D}$, and $\boldsymbol{D}_{d}$
        $\boldsymbol{v}_{K} \leftarrow \boldsymbol{R}_{K} \boldsymbol{v}$;
        $\boldsymbol{u}_{K} \leftarrow \boldsymbol{G}^{T}\left[\boldsymbol{D}_{d}\left(\boldsymbol{G} \boldsymbol{v}_{K}\right)\right]-k^{2} \boldsymbol{B}^{T}\left[\boldsymbol{D}\left(\boldsymbol{B} \boldsymbol{v}_{K}\right)\right] ;$
        $\boldsymbol{u} \leftarrow \boldsymbol{u}+\boldsymbol{R}_{K}^{T} \boldsymbol{u}_{K}$
```

to $O\left(P^{2 d}\right)$. Secondly, the product of a local matrix and vector must be computed in every iteration step [13], [14]. In Conjugate Gradient (CG) and GMRES methods, one SpMV multiplication is performed per iteration step, while twice are required for the BICGSTAB method [36].

## B. Basis Function Evaluation on Reference Cell

The Jacobian transformation at the $q$ th quadrature point $\boldsymbol{x}_{\boldsymbol{q}}$ is represented by a $d \times d$ Jacobian matrix

$$
\begin{equation*}
\boldsymbol{J}_{q}=\left.\frac{\partial \boldsymbol{x}}{\partial \boldsymbol{\xi}}\right|_{\boldsymbol{x}_{q}} \tag{16}
\end{equation*}
$$

where $\boldsymbol{\xi}=\left(\xi_{1}, \ldots, \xi_{d}\right)$ is the reference coordinate. To construct local matrices, one must initialize the basis function matrices and gradient matrices which are different on each cell. Instead of computing local matrices from (9) or (11), we evaluate the basis function values on the reference element [10], [25]

$$
\begin{align*}
& \left(\boldsymbol{W}_{K}\right)_{i j}=\sum_{q=1}^{\mathcal{Q}} \boldsymbol{J}_{q}^{-T} \nabla \hat{\Phi}_{i}\left(\boldsymbol{\xi}_{q}\right) \cdot \boldsymbol{J}_{q}^{-1} \nabla \hat{\Phi}_{j}\left(\boldsymbol{\xi}_{q}\right)\left|\boldsymbol{J}_{q}\right| w_{q} \\
& \left(\boldsymbol{M}_{K}\right)_{i j}=\sum_{q=1}^{\mathcal{Q}} \hat{\Phi}_{i}\left(\boldsymbol{\xi}_{q}\right) \hat{\Phi}_{j}\left(\boldsymbol{\xi}_{q}\right)\left|\boldsymbol{J}_{q}\right| w_{q} \tag{17}
\end{align*}
$$

or

$$
\begin{align*}
\boldsymbol{W}_{K} & =\hat{\boldsymbol{G}}^{T} \boldsymbol{J}_{K}^{-T} \hat{\boldsymbol{D}}_{d} \boldsymbol{J}_{K}^{-1} \hat{\boldsymbol{G}} \\
\boldsymbol{J}_{K} & =\left[\begin{array}{lll}
\boldsymbol{J}_{1}^{-1} & & \\
& \ddots & \\
& & \boldsymbol{J}_{\mathcal{Q}}^{-1}
\end{array}\right] \times \tilde{\boldsymbol{G}} \\
\boldsymbol{M}_{K} & =\hat{\boldsymbol{B}}^{T} \hat{\boldsymbol{D}} \hat{\boldsymbol{B}} \\
\hat{\boldsymbol{D}} & =\operatorname{diag}\left(w_{1}\left|\boldsymbol{J}_{1}\right|, w_{2}\left|\boldsymbol{J}_{2}\right|, \ldots, w_{\mathcal{Q}}\left|\boldsymbol{J}_{\mathcal{Q}}\right|\right) \\
\hat{\boldsymbol{D}}_{d} & =\hat{\boldsymbol{D}} \otimes \boldsymbol{I}_{d} \tag{18}
\end{align*}
$$

where $\hat{\Phi}, \hat{\boldsymbol{B}}, \hat{\boldsymbol{G}}, \hat{\boldsymbol{D}}$ and $\hat{\boldsymbol{D}}_{d}$ are corresponding basis function and matrices on reference cell. Their evaluation and storage are trivial since all physical elements with the same polynomial order share one reference cell. Therefore, only the determinant of the Jacobian matrix $\left|\boldsymbol{J}_{q}\right|$ and the inverse Jacobian matrix $\boldsymbol{J}_{q}^{-1}$ will be stored and the storage complexity is reduced from $O\left(P^{2 d}\right)$ to $O\left(P^{d}\right)$ [10].

## C. Sum Factorization for Lower MV Complexity

On the reference cell, the value and gradient of higher dimensional basis functions can be decomposed as the product

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-05.jpg?height=895&width=841&top_left_y=176&top_left_x=187)
Fig. 3. Decomposition of the second order basis functions on a reference cell. (a) Three-dimensional basis function $i$-indices are rearranged as (b) ( $i_{1}, i_{2}, i_{3}$ )-indices with second floor emitted. The value of a 3-D basis function is equivalent to the product of (c) 1-D basis functions.

of 1-D basis functions. For example, in hexahedral elements, the value and gradient of 3-D basis functions are

$$
\begin{align*}
& \hat{\Phi}_{i}\left(\boldsymbol{\xi}_{q}\right)=\hat{\Phi}_{i_{1} i_{2} i_{3}}\left(\boldsymbol{\xi}_{q}\right)=\phi_{i_{1}}\left(\xi_{1}^{(l)}\right) \phi_{i_{2}}\left(\xi_{2}^{(m)}\right) \phi_{i_{3}}\left(\xi_{3}^{(n)}\right) \\
& \left.\frac{\partial \hat{\Phi}_{i}}{\partial \xi_{1}}\right|_{\boldsymbol{\xi}_{q}}=\left.\frac{\partial \hat{\Phi}_{i_{1} i_{2} i_{3}}}{\partial \xi_{1}}\right|_{\boldsymbol{\xi}_{q}}=\phi_{i_{1}}^{\prime}\left(\xi_{1}^{(l)}\right) \phi_{i_{2}}\left(\xi_{2}^{(m)}\right) \phi_{i_{3}}\left(\xi_{3}^{(n)}\right) \tag{19}
\end{align*}
$$

where $\phi_{i_{j}}$ is the $i_{j}$ th 1-D basis function on $j$-direction [37, eq. (A10)], see Fig. 3. The unknown and quadrature indices are rearranged as $i=\left(i_{1}, i_{2}, i_{3}\right)$ and $\boldsymbol{\xi}_{\boldsymbol{q}}=\left(\xi_{1}^{(l)}, \xi_{2}^{(m)}, \xi_{3}^{(n)}\right)$, respectively. $l, m, n, i_{j} \in[1, P]$ since exact integration are performed. (19) can also be presented in an equivalent matrix form. The basis function matrix on the reference cell may be decomposed as the Kronecker product of 1-D basis function matrices

$$
\begin{equation*}
\hat{\boldsymbol{B}}=\boldsymbol{B}_{\xi_{1}} \otimes \boldsymbol{B}_{\xi_{2}} \otimes \boldsymbol{B}_{\xi_{3}} \tag{20}
\end{equation*}
$$

with

$$
\begin{equation*}
\left(\boldsymbol{B}_{\xi_{j}}\right)_{k i_{j}}=\phi_{i_{j}}\left(\xi_{j}^{(k)}\right), \quad k=l, m, \text { or } n . \tag{21}
\end{equation*}
$$

The same decomposition procedure is also applied to gradient matrices. Taking advantage of the tensor-product property of hexahedral elements, the sum factorization technique is applied to reduce the complexity of both local MV multiplication and local matrices assembly.

1) Local MV Evaluation With Sum Factorization: Generally, a local matrix-vector (MV) multiplication requires $O\left(P^{2 d}\right)$ complexity. On tensor-product elements, this process can be accelerated by the sum factorization technique. For example, the product of a mass matrix and the input vector in 2-D is equivalent to

$$
\hat{\boldsymbol{B}} \boldsymbol{v}_{K}=\left(\boldsymbol{B}_{\xi_{1}} \otimes \boldsymbol{B}_{\xi_{2}}\right) \boldsymbol{v}_{K}
$$

$$
\begin{equation*}
=\operatorname{vec}\left(\boldsymbol{B}_{\xi_{2}} \operatorname{vec}^{-1}\left(\boldsymbol{v}_{K}\right) \boldsymbol{B}_{\xi_{1}}^{T}\right) . \tag{22}
\end{equation*}
$$

where vec is the vectorization operator [38]. Equation (22) only has a complexity of $O\left(P^{d+1}\right)$. The multiplication between a Kronecker product and a vector in 3-D is expedited by rewriting the vector as a rank-3 tensor, see [35]. Because partial derivatives of basis functions also have tensor-product property, the strategy is also applied to other operators such as $(\nabla u, \nabla u)_{\Omega}$ and $(\nabla \times \boldsymbol{E}, \nabla \times \boldsymbol{E})_{\Omega}$. For scalar problems, a matrix-free global MV evaluation is faster than a matrix-assembled one when polynomial order is higher than 1 [10]. However, (4) is a vector problem and more basis functions are coupled in one element. An acceleration is observed in the matrix-free method even with linear basis functions, as shown in the numerical results.

The naming convention of sum factorization is more intuitive in non-Kronecker form. It involves decomposing the contributions of the basis functions to 1-D products and storing the intermediate value. For example, evaluate one mass matrix's entry as

$$
\begin{align*}
\left(\boldsymbol{M}_{K}\right)_{i j}= & \sum_{q=1}^{P^{3}} \hat{\Phi}_{i}\left(\boldsymbol{\xi}_{q}\right) \hat{\Phi}_{j}\left(\boldsymbol{\xi}_{q}\right) w_{q}\left|\boldsymbol{J}_{q}\right| \\
= & \sum_{l=1}^{P} \phi_{i_{1}}\left(\xi_{1}^{(l)}\right) \phi_{j_{1}}\left(\xi_{1}^{(l)}\right) w_{l} \\
& \times\left\{\sum_{m=1}^{P} \phi_{i_{2}}\left(\xi_{2}^{(m)}\right) \phi_{j_{2}}\left(\xi_{2}^{(m)}\right) w_{m}\right. \\
& \left.\quad\left[\sum_{n=1}^{P} \phi_{i_{3}}\left(\xi_{3}^{(n)}\right) \phi_{j_{3}}\left(\xi_{3}^{(n)}\right) w_{n}\left|\boldsymbol{J}_{l m n}\right|\right]\right\} . \tag{23}
\end{align*}
$$

Sum factorization does not accelerate the evaluation of individual entries, but simplifies the entire process of matrix assembly or local MV evaluation within one element. See [22] for detailed pseudocode of sum factorized local MV evaluation and further improvement leveraging the symmetric structure of $\boldsymbol{A}_{K}$.
2) Local Matrix-Assembly With Sum Factorization: An entry of $\boldsymbol{A}_{K}$ can be evaluated by a MV multiplication as

$$
\begin{equation*}
\left(\boldsymbol{A}_{K}\right)_{i j}=\left(\boldsymbol{A}_{K} \boldsymbol{e}^{j}\right)_{i} \tag{24}
\end{equation*}
$$

where $\boldsymbol{e}^{j}$ is the unit vector. Similarly, a column of $\boldsymbol{A}_{K}$ is computed as $\boldsymbol{A}_{K} \boldsymbol{e}^{j}$ [12]. Compared to traditional local FEM matrix assembly with $O\left(P^{3 d}\right)$ complexity, assembling $\boldsymbol{A}_{K}$ with sum factorization technique reduces the complexity to $O\left(P^{2 d}\right)$, as listed in Table I. Though $A_{K}$ will neither be assembled nor stored explicitly in the proposed algorithm, $A_{K}$ 's diagonal entries must be computed to run Jacobi smoother and multigrid preconditioner (detailed in Section IV). This strategy is used in many code projects [12], [13], [15], [39].

## IV. Matrix-Free Preconditioner Based on Multigrid Method

The linear system obtained from (6) has both positive and negative eigenvalues, therefore being an indefinite system [40]. It consists of both nonphysical negative eigenvalues

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-06.jpg?height=716&width=898&top_left_y=176&top_left_x=162)
Fig. 4. Global refined geometric multigrids on part of a 90 GHz integrated filter. (a) Geometry model. (b) Initial mesh. The geometry-embedded mesh refinement on (c) second level and (d) third level.

from irrotational space, and positive or positive and negative physical eigenvalues related to resonance frequency [41]. When $k_{0}$ coincides exactly with a resonant frequency, the system is singular [41]. Krylov subspace methods, such as CG and GMRES, do not guarantee convergence for indefinite or singular systems. Effective preconditioners are crucial to improve the convergence of Krylov subspace methods. The main difficulty is the matrix-free scheme limits our option when selecting a preconditioner, since popular preconditioners like successive over-relaxation (SOR), incomplete LU (ILU) decomposition, and algebraic multigrid (AMG) require $\boldsymbol{A}$ be stored explicitly.

A potential candidate is the Jacobi preconditioner because the diagonal entries are cheap to compute with element-wise MV multiplications as shown in (24). The global diagonal is constructed from element-wise diagonals with $\boldsymbol{R}_{K}^{T}$ operators. However, Jacobian preconditioned systems need more iterations than SOR and ILU [36]. Another preconditioner available in the matrix-free scheme is the geometric multigrid preconditioner [17], [20], [21]. Multigrid or multilevel preconditioners suppress both high-frequency and low-frequency residuals. A direct solver is applied on the sufficiently coarse level. The polynomial multigrid is a popular option for commercial software such as Ansys HFSS [42] and COMSOL Multiphysics ${ }^{1}$ [43]. In a matrix-free framework, the matrix on the coarsest level will be assembled explicitly.

In this work, we implement a matrix-free geometric multigrid preconditioner with Chebyshev iteration smoothers [44]. The diagonals on each level are computed with a forward matrix-free evaluation as (24). GMRES with a restart number of 50 is chosen as the outer Krylov subspace method for all numerical examples. We choose an octree-like refinement strategy where one macro element is divided into $2^{d}$ subelements, as shown in Fig. 4. Geometry information is used to

[^1]![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-06.jpg?height=586&width=776&top_left_y=184&top_left_x=1130)
Fig. 5. Scattering analysis of a PEC sphere with 0.2 m radius. The real part of electric field in $z$-direction is shown.

guide curvilinear multigrid refinement and reduce discretization error when combined with high-order basis functions [45].

Despite the capability of the current solver configuration to yield satisfactory results as will be demonstrated in Section V, more appropriate preconditioners are recommended for indefinite Helmholtz equation or time-harmonic Maxwell's wave equation. Otherwise, the coarsest level must be fine enough to satisfy the points per wavelengths (PPWs) constraint for wave problems [46]. For higher wave numbers, even constant $k h$ product is not sufficient to avoid the pollution effect, and relative error in $H^{1}$-norm is governed by $k^{2} h^{3}$ [47]. Common multilevel preconditioners include domain decomposition method [8], [48], shifted laplacian preconditioner [49], [50] and sweeping preconditioner [51], [52], etc. Their implementation under a matrix-free framework necessitates further investigation. Some matrix-free auxiliary space preconditioners have been proposed for positive definite Maxwell system [53], [54], [55].

## V. Numerical Results

In this section, a sequence of examples is executed to showcase the reduced memory and runtime properties of the proposed matrix-free electromagnetic solver, relative to the conventional matrix-assembled solver. All parallel computation is performed with distributed memory. P4est is an open-source library for adaptive mesh refinement on forests of octrees [56]. Our work utilizes p4est for the parallel mesh partition. All examples are tested on a workstation with 4 Intel Xeon Platinum 8280M CPUs at 2.70 GHz and 1.2 TB memory.

## A. Wave Scattering From a Metal Sphere

This example validates the correctness of our matrix-free FEM solver on a metal sphere with a $0.2-\mathrm{m}$ radius, radiated by a 2 GHz z-polarized plane wave propagating from the $-x$ direction to the $+x$-direction. In this example, a two-level multigrid preconditioner and second-order basis functions are applied. Fig. 5 displays the electric fields in the $z$ direction. Bistatic radar cross sections (RCSs) computed by the matrix-free method and analytical MIE scattering method

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-07.jpg?height=595&width=796&top_left_y=178&top_left_x=205)
Fig. 6. RCS of the metal sphere at 2 GHz computed by MIE approach and matrix-free FEM, respectively. The relative error between numerical and analytical solutions is $1.4 \%$.

are shown in Fig. 6. The analytical and numerical RCS results are in good agreement in the strong reflection area, and the difference is less than 1.1 dB in the eclipsed area. To reduce discretization error, the second-order curvilinear mesh is applied to capture the spherical surface.

## B. Analytical Waveguide

This example verifies the mesh convergence of the proposed method. Matrix-free FEM is applied on a $0.03 \times 0.015 \times$ 0.45 m rectangular waveguide to simulate its $\mathrm{TE}_{20}$ electric field as shown in Fig. 7(a). The accuracy of FEM approximation is verified by the relative $L^{2}$ error

$$
\begin{equation*}
e=\frac{\|\boldsymbol{E}(\boldsymbol{r})-\hat{\boldsymbol{E}}(\boldsymbol{r})\|_{L^{2}}}{\|\hat{\boldsymbol{E}}(\boldsymbol{r})\|_{L^{2}}} \tag{25}
\end{equation*}
$$

where $\hat{\boldsymbol{E}}(\boldsymbol{r})$ is the analytical electric field. Fig. 7(b) gives the exponential mesh convergence for up to third polynomial order. The optimal $p+1$ convergence is achieved for Lagrange basis functions, which demonstrates the reliability of the proposed solver.

## C. 16-Way High-Power Waveguide Divider

This example demonstrates the superior efficiency and scalability of our algorithm on a 16-way waveguide power divider as compared to conventional FEM. These structures are typically used for power transfer from high-power microwave sources to antenna arrays. Electric fields in the $z$-direction are depicted in Fig. 8 where halved power is observed at each junction. Fig. 9 presents the S -parameters from 15 to 17 GHz computed by the proposed algorithm and COMSOL Multiphysics ${ }^{1}$ [43]. A good agreement between S-parameters is achieved. Since S-parameters are more sensitive criteria than electric fields, it is sufficient to prove the accuracy of the proposed matrix-free algorithm. The in-house solver is paralleled through OpenMP with fully distributed memory, while COMSOL supports hybrid (MPI and OpenMP) parallelism.

Strong scalability tests at 16 GHz are performed for a four-level multigrid case (around 43 million DoF ) and a

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-07.jpg?height=747&width=858&top_left_y=178&top_left_x=1086)
Fig. 7. Mesh convergence test for the analytical waveguide. (a) Real part of electric field in the rectangular waveguide and (b) exponential drops of numerical error. Optimal mesh convergence is achieved.

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-07.jpg?height=266&width=866&top_left_y=1073&top_left_x=1081)
Fig. 8. Electric field distribution in the $z$-direction of the 16 -way power divider.

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-07.jpg?height=627&width=871&top_left_y=1466&top_left_x=1086)
Fig. 9. S-parameters of the 16-way power divider from 15 to 17 GHz .

five-level multigrid case (around 334 million DoF). Secondorder basis functions and tetrahedral mesh are used in COMSOL. Fig. 10 provides the solution time of different methods. For four-level multigrid, the mesh size is $1 / 24$ wavelength, and the memory cost is 22,520 , and 180 GB for the matrix-free solver, in-house matrix-assembled solver, and COMSOL, respectively. COMSOL's frequency-domain RF solver is faster than our in-house matrix-assembled solver with no more than 48 cores. However, it seldom gains acceleration

TABLE II
Memory and Runtime Statistics of Power Divider by Matrix-Free (MF) and Matrix-Assembled (MA) Methods
| Algo. \& $p$ | Mesh (GB) | $\boldsymbol{A}$ (GB) | $\boldsymbol{A}_{M}$ (GB) | $\boldsymbol{J}_{q}^{-1}$ (GB) | Total (GB) | Finest level assembly (s) | Multigrid assembly (s) | Iteration (s) | GMRES Steps | Total (s) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| MA \& 2 | 0.156 | 24.1 | 27.61 | - | 78.2 | 21.5 | 25.5 | 57.4 | 36 | 110.0 |
| MF \& 2 | 0.156 | - | 0.4 | 0.61 | 9.04 | 0 | 0.5 | 13.8 | 34 | 14.6 |
| MA \& 1 | 0.565 | 10.3 | 11.8 | - | 46.4 | 5.0 | 5.9 | 42.3 | 39 | 52.9 |
| MF \& 1 | 0.565 | - | 0.178 | 1.38 | 9.99 | 0 | 0.315 | 23.8 | 38 | 25.5 |


$\boldsymbol{A}, \boldsymbol{A}_{M}$ and $\boldsymbol{J}_{q}^{-1}$ refers to the global matrix, multigrid matrices, and the inverse Jacobian matrices, respectively. The multigrid matrices storage of matrix-free methods is not zero because the matrix on the coarsest level is assembled.

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-08.jpg?height=638&width=827&top_left_y=672&top_left_x=195)
Fig. 10. Strong scalability test for the 16 -way waveguide. The mesh is refined to four levels (43 million DoF, square labeled) and five levels ( 334 million DoF, triangle labeled), respectively. Matrix-free solver consumes less time and has better scalability than COMSOL.

from increased CPU resources because COMSOL applies a SOR smoother. SOR preconditioner is essentially a sequential algorithm whose parallelism decreases with more sub-block numbers [36], [44]. Compared to the matrix-assembled methods (both in-house and COMSOL), the matrix-free method exhibits a significant reduction in memory and time consumption. The speedup decreases as the computation time reaches 60 s or less.

When the multigrid is refined to five levels, the limited memory resources prohibit the test of traditional FEM. However, the matrix-free method is still available and exhibits an even shorter runtime than matrix-assembled solver with fourlevel multigrids. Our proposed algorithm exhibits excellent strong scalability with less than or equal to 48 cores.

We further compare the detailed time and memory consumption for both first and second-order polynomial basis functions. Table II provides the result summarization run with 40 cores. To keep the number of DoFs the same, the mesh level is 4 and 3 for first and second-order polynomials, respectively. For the second-order polynomials, matrix assembly takes up more than 65\% memory and 40\% time of the FEM solver with matrixassembled method. By comparison, with matrix-free method, the GMRES solver takes up $94 \%$ of the total runtime. Another distinctive advantage of the matrix-free method is that the higher polynomial order requires less computation time than

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-08.jpg?height=340&width=890&top_left_y=677&top_left_x=1073)
Fig. 11. Cross section of the tex-hex mesh of the missile. The number of DoF in this problem amounts to 22.5 million.

the lower one, as opposed to the matrix-assembled FEM. With an equal number of unknowns, the computation reduction from cell numbers outweighs the increment from the basis function order. Thus, the proposed algorithm is highly applicable to high-order FEMs. Finally, the GMRES solving stage is also accelerated by the single instruction multiple data (SIMD) vectorization technique [10], [57]. It allows a double array consisting of eight doubles to be processed in one instruction because our CPU supports the AVX512 instruction set. This technique is hardly realizable for SpMV multiplication in conventional FEMs.

## D. Wave Scattering From a Missile

The final example involves the electric field scattered from a missile. Meshes consisting of hexahedra elements are difficult to generate for complex geometries. As a meshing strategy, we first divide the domain into tetrahedra and then partition each tetrahedron into 4 hexahedra. A post-mesh-smoothing process is deemed necessary due to the low mesh quality for partitioned hexahedra and the possibility of negative element Jacobian around surface corners. The negative impact of poor mesh quality on time-domain methods should not be overlooked, as it induces instability. However, the adverse effect on frequency-domain problems is relatively slight, which only increases the conditional number and entails additional Krylov subspace iterations.

The missile is 4 m long with metal material coverings, illuminated by a y-polarized plane wave from nose to tail at 2.4 GHz . We employ third-order polynomials and second-order curvilinear mesh on a two-level multigrid, resulting in 22.5 million freedoms. Computation is reduced by leveraging the geometric symmetry through PEC and perfect magnetic conductor (PMC) conditions.

![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-09.jpg?height=906&width=895&top_left_y=181&top_left_x=157)
Fig. 12. Scattering analysis of a missile-like object. The length of the body is 4 m . A 2.4 GHz plane wave polarized in $y$-direction is incident from $+z$ to $-z$. (a) Real part of the electric field in the $y$-direction and (b) RCS computed by matrix-free method and FEKO.

TABLE III
Runtime and Memory of Missile Scattering Problem With Matrix-Assembled Method and Matrix-Free Method
|  | Memory | Assembly time | Total time |
| :---: | :---: | :---: | :---: |
| MA | 947 GB | 764 s | 2913 s |
| MF | 439 GB | - | 1394 s |


In this work, the LU decomposition is realized by MUMPS solver [58]. The electric field on the xOz cross section is depicted in Fig. 12(a) and 12(b) compares the RCS results computed by the matrix-free solver and FEKO [59]. Table III demonstrates the total time and memory consumption of both methods. The improvement from matrix-free method is not as apparent as in Table II due to the expensive LU decomposition of a large coarse-level matrix. It consumes around 650 s runtime and 350 GB of memory. Performance improvement is decreased compared to the waveguide example because the computation time and memory are dominated by the direct solver. However, the matrix-free method still has a significant improvement over the traditional FEM.

## VI. Conclusion

This article presents a fast and low-storage matrix-free algorithm for solving the time-harmonic Maxwell's wave equation. By decomposing global SpMV multiplications at a cellular level, our approach circumvents the need for the assembly of a global matrix that is necessary for matrixassembled FEMs. Additionally, the use of sum factorization technique on tensorial basis functions assists in reducing the computational cost associated with local MV evaluation.

A geometric multigrid preconditioner is applied to improve the convergence of the indefinite linear system discretized from the time-harmonic Maxwell's equations. Numerical examples demonstrate its fast and low-storage properties and good scalability. The proposed method is also friendly to high-order methods due to reduced computational and storage complexities w.r.t. the polynomial order of basis functions. Future work includes extending the matrix-free frequency-domain FEM to Nédélec elements and efficient EM-thermal-mechanical multiphysics coupling numerical methods.

## References

[1] Q. Zhan, M. Zhuang, Y. Mao, and Q. H. Liu, "Unified Riemann solution for multi-physics coupling: Anisotropic poroelastic/elastic/fluid interfaces," J. Comput. Phys., vol. 402, Feb. 2020, Art. no. 108961.
[2] Z. Peng, K.-H. Lim, and J.-F. Lee, "Nonconformal domain decomposition methods for solving large multiscale electromagnetic scattering problems," Proc. IEEE, vol. 101, no. 2, pp. 298-319, Feb. 2013.
[3] J. Wang and W.-Y. Yin, "Development of a novel FDTD (2, 4)compatible conformal scheme for electromagnetic computations of complex curved PEC objects," IEEE Trans. Antennas Propag., vol. 61, no. 1, pp. 299-309, Jan. 2013.
[4] J.-M. Jin, The Finite Element Method in Electromagnetics. Hoboken, NJ, USA: Wiley, 2015.
[5] Q. Zhan, Y. Fang, M. Zhuang, M. Yuan, and Q. H. Liu, "Stabilized DGPSTD method with nonconformal meshes for electromagnetic waves," IEEE Trans. Antennas Propag., vol. 68, no. 6, pp. 4714-4726, Jun. 2020.
[6] Q. Zhan et al., "An adaptive high-order transient algorithm to solve largescale anisotropic Maxwell's equations," IEEE Trans. Antennas Propag., vol. 70, no. 3, pp. 2082-2092, Mar. 2022.
[7] P.-H. Jia et al., "Twofold domain decomposition method for the analysis of multiscale composite structures," IEEE Trans. Antennas Propag., vol. 67, no. 9, pp. 6090-6103, Sep. 2019.
[8] H.-X. Zhang et al., "Multiphysics computing of challenging antenna arrays under a supercomputer framework," IEEE J. Multiscale Multiphys. Comput. Tech., vol. 8, pp. 165-177, 2023.
[9] S. A. Orszag, "Spectral methods for problems in complex geometrics," in Numer. Methods Partial Differ. Equ. Amsterdam, The Netherlands: Elsevier, 1979, pp. 273-305.
[10] M. Kronbichler and K. Kormann, "A generic interface for parallel cellbased finite element operator application," Comput. Fluids, vol. 63, pp. 135-147, Jun. 2012.
[11] G. Schubert, G. Hager, and H. Fehske, "Performance limitations for sparse matrix-vector multiplications on current multi-core environments," in Proc. High Perform. Comput. Sci. Eng. Trans. 4th Joint HLRB KONWIHR Rev. Results Workshop, Leibniz Supercomputing Centre, Garching/Munich, Germany. Springer, Dec. 2009, pp. 13-26.
[12] M. Kronbichler and K. Ljungkvist, "Multigrid for matrix-free highorder finite element computations on graphics processors," ACM Trans. Parallel Comput., vol. 6, no. 1, pp. 1-32, Mar. 2019.
[13] D. Arndt et al., "ExaDG: High-order discontinuous Galerkin for the exa-scale," in Proc. Softw. Exascale Comput.-SPPEXA, H.-J. Bungartz, S. Reiz, B. Uekermann, P. Neumann, and W. E. Nagel, Eds. Springer Int. Publishing, 2016, pp. 189-224.
[14] P. Fischer et al., "NekRS, a GPU-accelerated spectral element NavierStokes solver," Parallel Comput., vol. 114, Dec. 2022, Art. no. 102982.
$[15]$ D. Arndt et al., "The deal.II library, Version 9.4," J. Numer. Math., vol. 30, no. 3, pp. 231-246, 2022.
[16] F. Dugast et al., "Part-scale thermal process modeling for laser powder bed fusion with matrix-free method and GPU computing," Additive Manuf., vol. 37, Jan. 2021, Art. no. 101732.
[17] J. Rudi et al., "An extreme-scale implicit solver for complex PDEs: Highly heterogeneous flow in Earth's mantle," in Proc. Int. Conf. High Perform. Comput. Netw. Storage Anal., Nov. 2015, pp. 1-12.
[18] J. Martínez-Frutos, P. J. Martínez-Castejón, and D. Herrero-Pérez, "Finegrained GPU implementation of assembly-free iterative solver for finite element problems," Comput. Struct., vol. 157, pp. 9-18, Sep. 2015.
[19] M. Kronbichler and K. Kormann, "Fast matrix-free evaluation of discontinuous Galerkin finite element operators," ACM Trans. Math. Softw., vol. 45, no. 3, pp. 1-40, Sep. 2019.
[20] P. Munch, I. Dravins, M. Kronbichler, and M. Neytcheva, "Stage-parallel fully implicit Runge-Kutta implementations with optimal multilevel preconditioners at the scaling limit," 2022, arXiv:2209.06700.
$[21]$ P. Munch, T. Heister, L. Prieto Saavedra, and M. Kronbichler, "Efficient distributed matrix-free multigrid methods on locally refined meshes for FEM computations," ACM Trans. Parallel Comput., vol. 10, no. 1, pp. 1-38, Mar. 2023.
[22] J. Mora and L. Demkowicz, "Fast integration of DPG matrices based on sum factorization for all the energy spaces," Comput. Methods Appl. Math., vol. 19, no. 3, pp. 523-555, Jul. 2019.
[23] N. Krais et al., "FLEXI: A high order discontinuous Galerkin framework for hyperbolic-parabolic conservation laws," Comput. Math. Appl., vol. 81, pp. 186-219, Jan. 2021.
$[24]$ D. Moxey et al., "Nektar++: Enhancing the capability and application of high-fidelity spectral/hp element methods," Comput. Phys. Commun., vol. 249, Apr. 2020, Art. no. 107110.
[25] D. Moxey, R. Amici, and M. Kirby, "Efficient matrix-free high-order finite element evaluation for simplicial elements," SIAM J. Sci. Comput., vol. 42, no. 3, pp. C97-C123, Jan. 2020.
[26] T. Sun, L. Mitchell, K. Kulkarni, A. Klöckner, D. A. Ham, and P. H. Kelly, "A study of vectorization for matrix-free finite element methods," Int. J. High Perform. Comput. Appl., vol. 34, no. 6, pp. 629-644, Nov. 2020.
[27] Z. Peng and D. Appelö, "EM-WaveHoltz: A flexible frequency-domain method built from time-domain solvers," IEEE Trans. Antennas Propag., vol. 70, no. 7, pp. 5659-5671, Jul. 2022.
[28] M. S. Tong and W. C. Chew, "Multilevel fast multipole acceleration in the nyström discretization of surface electromagnetic integral equations for composite objects," IEEE Trans. Antennas Propag., vol. 58, no. 10, pp. 3411-3416, Jul. 2010.
[29] E. Agullo, L. Giraud, A. Gobé, M. Kuhn, S. Lanteri, and L. Moya, "High order HDG method and domain decomposition solvers for frequencydomain electromagnetics," Int. J. Numer. Modelling, Electron. Netw. Devices Fields, vol. 33, no. 2, p. e2678, Mar. 2020.
[30] H. Feng, T.-Y. Li, M. Zhuang, H. Xie, W.-Y. Yin, and Q. Zhan, "Skeleton-enhanced discontinuous Galerkin method for 3-D nonlinear semiconductor modeling," IEEE Trans. Microwave Theory Techn., vol. 71, no. 6, pp. 1-13, Jun. 2023.
[31] S. Petrides and L. F. Demkowicz, "An adaptive DPG method for high frequency time-harmonic wave propagation problems," Comput. Math. Appl., vol. 74, no. 8, pp. 1999-2017, Oct. 2017.
[32] W. H. Press, W. T. Vetterling, S. A. Teukolsky, and B. P. Flannery, Numer. Recipes. Princeton, NJ, USA: Citeseer, 1988.
[33] J.-C. Nédélec, "Mixed finite elements in $\mathbb{R}^{3}$," Numerische Math., vol. 35, pp. 315-341, 1980.
[34] P. E. J. Vos, S. J. Sherwin, and R. M. Kirby, "From h to p efficiently: Implementing finite and spectral/hp element methods to achieve optimal performance for low- and high-order discretisations," J. Comput. Phys., vol. 229, no. 13, pp. 5161-5181, Jul. 2010.
[35] C. D. Cantwell, S. J. Sherwin, R. M. Kirby, and P. H. J. Kelly, "From h to p efficiently: Strategy selection for operator evaluation on hexahedral and tetrahedral elements," Comput. Fluids, vol. 43, no. 1, pp. 23-28, Apr. 2011.
[36] Y. Saad, Iterative Methods for Sparse Linear Systems. New Delhi, India: SIAM, 2003.
[37] Q. Zhan, Q. Sun, Q. Ren, Y. Fang, H. Wang, and Q. H. Liu, "A discontinuous Galerkin method for simulating the effects of arbitrary discrete fractures on elastic wave propagation," Geophys. J. Int., vol. 210, no. 2, pp. 1219-1230, Aug. 2017.
[38] P. E. Buis and W. R. Dyksen, "Efficient vector and parallel manipulation of tensor products," ACM Trans. Math. Softw., vol. 22, no. 1, pp. 18-23, Mar. 1996.
[39] M. Homolya, R. C. Kirby, and D. A. Ham, "Exposing and exploiting structure: Optimal code generation for high-order finite element methods," 2017, arXiv:1711.02473.
[40] D.-K. Sun, J.-F. Lee, and Z. Cendes, "Construction of nearly orthogonal Nedelec bases for rapid convergence with multilevel preconditioned solvers," SIAM J. Sci. Comput., vol. 23, no. 4, pp. 1053-1076, Jan. 2001.
[41] R. Dyczij-Edlinger and O. Biro, "A joint vector and scalar potential formulation for driven high frequency problems using hybrid edge and nodal finite elements," IEEE Trans. Microwave Theory Techn., vol. 44, no. 1, pp. 15-23, Jan. 1996.
[42] Ansys HFSS Softw. Accessed: Feb. 8, 2024. [Online]. Available: https://www.ansys.com/products/electronics/ansys-hfss
[43] COMSOL Multiphysics V.6.2, COMSOL AB, Stockholm, Sweden. Accessed: Feb. 8, 2024. [Online]. Available: https://www.comsol.com
[44] M. Adams, M. Brezina, J. Hu, and R. Tuminaro, "Parallel multigrid smoothing: Polynomial versus Gauss-Seidel," J. Comput. Phys., vol. 188, no. 2, pp. 593-610, Jul. 2003.
[45] L. Heltai, W. Bangerth, M. Kronbichler, and A. Mola, "Propagating geometry information to finite element computations," ACM Trans. Math. Softw., vol. 47, no. 4, pp. 1-30, Dec. 2021.
[46] K. H. Dridi, J. S. Hesthaven, and A. Ditkowski, "Staircase-free finitedifference time-domain formulation for general materials in complex geometries," IEEE Trans. Antennas Propag., vol. 49, no. 5, pp. 749-756, May 2001.
[47] F. Ihlenburg and I. Babuška, "Finite element solution of the Helmholtz equation with high wave number Part I: The h-version of the FEM," Comput. Math. Appl., vol. 30, no. 9, pp. 9-37, Nov. 1995.
[48] K. Zhao, V. Rawat, and J.-F. Lee, "A domain decomposition method for electromagnetic radiation and scattering analysis of multi-target problems," IEEE Trans. Antennas Propag., vol. 56, no. 8, pp. 2211-2221, Aug. 2008.
[49] Y. A. Erlangga, C. Vuik, and C. W. Oosterlee, "On a class of preconditioners for solving the Helmholtz equation," Appl. Numer. Math., vol. 50, nos. 3-4, pp. 409-425, Sep. 2004.
[50] G. Rizzuti and W. A. Mulder, "Multigrid-based-shifted-Laplacian' preconditioning for the time-harmonic elastic wave equation," J. Comput. Phys., vol. 317, pp. 47-65, Jul. 2016.
[51] B. Engquist and L. Ying, "Sweeping preconditioner for the Helmholtz equation: Moving perfectly matched layers," Multiscale Model. Simul., vol. 9, no. 2, pp. 686-710, Apr. 2011.
[52] P. Tsuji, B. Engquist, and L. Ying, "A sweeping preconditioner for timeharmonic Maxwell's equations with finite elements," J. Comput. Phys., vol. 231, no. 9, pp. 3770-3783, May 2012.
[53] W. Pazner, "Efficient low-order refined preconditioners for high-order matrix-free continuous and discontinuous Galerkin methods," SIAM J. Sci. Comput., vol. 42, no. 5, pp. A3055-A3083, Jan. 2020.
[54] A. T. Barker and T. Kolev, "Matrix-free preconditioning for highorderH(curl) discretizations," Numer. Linear Algebra Appl., vol. 28, no. 2, p. e2348, Mar. 2021.
[55] W. Pazner, T. Kolev, and C. R. Dohrmann, "Low-order preconditioning for the high-order finite element de Rham complex," SIAM J. Sci. Comput., vol. 45, no. 2, pp. A675-A702, Apr. 2023.
[56] C. Burstedde, L. C. Wilcox, and O. Ghattas, "P4est: Scalable algorithms for parallel adaptive mesh refinement on forests of octrees," SIAM J. Sci. Comput., vol. 33, no. 3, pp. 1103-1133, Jan. 2011.
[57] B. Bramas and P. Kus, "Computing the sparse matrix vector product using block-based kernels without zero padding on processors with AVX-512 instructions," PeerJ Comput. Sci., vol. 4, p. e151, Apr. 2018.
[58] P. Amestoy, I. S. Duff, J.-Y. L'Excellent, and J. Koster, "MUltifrontal massively parallel solver (MUMPS Version 4.2 beta) user's guide," Ph.D. dissertation, Dept. Laboratoire de l'informatique du parallélisme, 2002.
[59] FEKO. Accessed: Feb. 8, 2024. [Online]. Available: https://altair. com/feko
![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-10.jpg?height=329&width=266&top_left_y=2276&top_left_x=1076)

Yi-Yao Wang (Graduate Student Member, IEEE) received the B.E. degree in information engineering from the College of Information Science and Electronic Engineering, Zhejiang University, Hangzhou, China, in 2021, where he is currently pursuing the Eng.D. degree in electronic science and technology.

His research interests include finite element method, multiphysics modeling, high-performance electromagnetic computing, and parallel computation.
![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-11.jpg?height=327&width=264&top_left_y=205&top_left_x=165)

Qiwei Zhan (Member, IEEE) received the B.S. degree from the University of Science and Technology of China, Hefei, China, in 2013, and the M.S. degree in civil and environmental engineering and the Ph.D. degree in electrical and computer engineering from Duke University, Durham, NC, USA, in 2016 and 2019, respectively.

From June 2019 to August 2020, he was a Peter O'Donnell, Jr. Post-Doctoral Fellow at the Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin, Austin, TX, USA. Since September 2020, he has been with the College of Information Science and Electronic Engineering, Zhejiang University, Hangzhou, China, as a tenure-track Professor and a Ph.D. Supervisor. His research interests include multiphysics modeling, computational mechanical waves, computational electromagnetics, uncertainty quantification, effective medium theory, inverse problems, biomedical imaging, and scientific machine learning.
![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-11.jpg?height=324&width=264&top_left_y=861&top_left_x=162)

Haoqiang Feng was born in Zhejiang, China, in 1997. He received the B.E. degree in information engineering from Zhejiang University, Hangzhou, China, in 2019, where he is currently pursuing the Ph.D. degree in electronic science and technology.

His research interests include advanced numerical methods, multiphysics modeling and highperformance computing.
![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-11.jpg?height=321&width=261&top_left_y=1303&top_left_x=165)

Yin-Da Wang (Student Member, IEEE) received the B.E. degree in information engineering from the Nanjing University of Aeronautics and Astronautics, Nanjing, China, in 2020. He is currently pursuing the Ph.D. degree in electronic science and technology with Zhejiang University, Hangzhou, China.

His research interests include finite element method, multiphysics modeling, computational electromagnetics, optimization theory, and highperformance computing.
![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-11.jpg?height=327&width=263&top_left_y=205&top_left_x=1076)

Bozhao Sun received the B.S. degree from Zhejiang University, Hangzhou, China, where he is currently pursing the Ph.D. degree.

His research interests include uncertainty quantification, inverse design, and scientific machine learning.
![](https://cdn.mathpix.com/cropped/2025_09_10_ca73905fcc086da7bde3g-11.jpg?height=328&width=273&top_left_y=865&top_left_x=1072)

Wen-Yan Yin (Fellow, IEEE) received the M.Sc. degree in electromagnetic field and microwave technique from Xidian University, Xi'an, China, in 1989, and the Ph.D. degree in electrical engineering from Xi'an Jiaotong University, Xi'an, in 1994.

He was an Associate Professor with the Department of Electronic Engineering, Northwestern Polytechnic University, Xi'an, from 1993 to 1996. From 1996 to 1998, he was a Research Fellow with the Department of Electrical Engineering, Duisburg University, Duisburg, Germany. Since 1998, he has been with the Monolithic Microwave Integrated Circuit Modeling and Package Laboratory, Department of Electrical Engineering, National University of Singapore (NUS), Singapore, as a Research Fellow. In 2002, he joined Temasek Laboratories, NUS, as a Research Scientist and the Project Leader of High-Power Microwave and Ultrawideband Electromagnetic Compatibility (EMC)/Electromagnetic Interference. Since 2005, he has been a Professor with the School of Electronic Information and Electrical Engineering, Shanghai Jiao Tong University, Shanghai, China. In 2009, he joined the National Key Laboratory of Modern Optical Instrumentation, Zhejiang University, Hangzhou, China, as a Qiu Shi Distinguished Professor. His current research interests include passive and active RF and millimeter-wave device and circuit modeling, ultrawideband interconnects and signal integrity, nanoelectronics, EMC and electromagnetic protection of communication platforms, and computational multiphysics and its application.


[^0]:    Manuscript received 3 July 2023; revised 9 January 2024; accepted 24 January 2024. Date of publication 6 February 2024; date of current version 7 March 2024. This work was supported in part by the National Key Research and Development Program of China under Grant 2023YFB3905003; in part by the National Natural Science Foundation of China under Grant 62071418, Grant 61931007, and Grant 92066105; and in part by the Zhejiang Provincial Natural Science Foundation of China under Grant LDT23F04014F01. (Corresponding authors: Qiwei Zhan; Wen-Yan Yin.)

    The authors are with the Innovative Institute of Electromagnetic Information and Electronic Integration, Key Laboratory of Advanced Micro-Nano Electronic Devices and Smart Systems, Zhejiang University, Hangzhou 310027, China (e-mail: qwzhan@zju.edu.cn; wyyin@zju.edu.cn).

    Color versions of one or more figures in this article are available at https://doi.org/10.1109/TAP.2024.3360700.

    Digital Object Identifier 10.1109/TAP.2024.3360700

[^1]:    ${ }^{1}$ Registered trademark.

