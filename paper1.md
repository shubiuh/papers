# Critical Review: "Differentiable finite element method with Galerkin discretization for fast and accurate inverse analysis of multidimensional heterogeneous engineering structures"

## Paper Summary

**Authors**: Xi Wang, Zhen-Yu Yin, Wei Wu, He-Hua Zhu  
**Journal**: Computer Methods in Applied Mechanics and Engineering (2025)  
**DOI**: [10.1016/j.cma.2025.117755](https://doi.org/10.1016/j.cma.2025.117755)

### Main Claims
- Proposes a "novel" Differentiable Finite Element Method (DFEM) 
- Claims to be a Physics-Encoded Numerical Network (PENN)
- Reports 3-5 orders of magnitude improvement over PINNs
- Demonstrates inverse analysis for 3D engineering structures

## Critical Analysis

### ❌ **Novelty Claims Are Overstated**

#### What the Paper Claims:
> "Novel differentiable finite element method based on Galerkin discretization"
> "Physics-Encoded Numerical Network (PENN)"

#### Reality Check:
**FEM + Automatic Differentiation has existed since the 1980s:**

- **Dixon et al. (1988)**: "Finite element optimization: The use of structured automatic differentiation"
- **Multiple papers throughout 1990s-2010s** on AD for FEM sensitivity analysis
- **MOOSE framework** integrated AD in 2018
- **deal.II, MFEM** have had AD capabilities for years

**What DFEM actually is:**
```
DFEM = Standard FEM + Automatic Differentiation + PyTorch/TensorFlow
```

No neural networks involved - just gradient-based parameter optimization.

### ⚠️ **Misleading Terminology**

#### "Differentiable Computational Graph"
This is just standard automatic differentiation that any modern optimization framework provides. Every time you call `loss.backward()` in PyTorch, you're using a "differentiable computational graph."

#### "Physics-Encoded Numerical Network" 
There's no neural network architecture. It's classical finite element discretization with AD-computed gradients.

#### "Neural Network Training"
The only "training" is parameter optimization: `θ ← θ - α∇L(θ)` where θ are physical parameters (material properties, boundary conditions).

### 📏 **Scale of Examples is Insufficient**

| Example | Parameters | Elements | Complexity |
|---------|------------|----------|------------|
| Thick cylinder | 1 (pressure) | ~1000 | Toy problem |
| Indentation | 2 (material params) | ~2000 | Toy problem |
| Cantilever beam | 20-40 (materials) | ~10000 | Still toy problem |
| Tunnel lining | 2 (loadings) | ~1000 | Toy problem |

**Industry-scale problems typically involve:**
- **10³ - 10⁶ parameters** (e.g., distributed material properties)
- **10⁶ - 10⁹ elements** (reservoir models, large structures)
- **Real-time constraints** (seconds to minutes)

### 📈 **Performance Issues with Realistic Data**

#### Convergence with Partial Observations
The tunnel lining example (Section 5.3) shows concerning behavior:
- **1,079 iterations** to converge
- **13 seconds** for a toy problem with only 2 parameters
- Only 1/3 of nodes had measurement data

**Scaling implications:**
- For 1000 parameters: ~500,000 iterations?
- For 10⁶ parameters: Computationally intractable

#### Memory and Computational Scaling
- **AD memory overhead**: 2-5× forward simulation
- **Gradient computation**: O(n_params) for forward-mode AD
- **No discussion** of reverse-mode AD or computational complexity

### 🔍 **Missing Critical Comparisons**

The paper only compares against PINNs, but ignores established methods for inverse problems:

#### Industry Standard Methods:
- **Ensemble Kalman Filter (EnKF)** - handles uncertainty, no gradients needed
- **Adjoint methods** - classical approach for large-scale inverse problems  
- **Reduced-order models** - enable real-time inversion
- **Gradient-free optimization** - genetic algorithms, particle swarms
- **Machine learning surrogates** - pre-trained fast approximations

#### Why These Matter:
- Handle **uncertainty quantification**
- Scale to **realistic problem sizes**
- Provide **real-time performance**
- **Don't require expensive gradient computations**

### 🎯 **What the Paper Actually Contributes**

#### Positive Aspects:
1. **Clean implementation** of FEM+AD in modern ML frameworks
2. **Good benchmarking** against PINNs (shows classical methods can outperform)
3. **Demonstrates value** of physics-informed approaches vs pure ML
4. **Educational value** - shows how to implement FEM in PyTorch/TensorFlow

#### Actual Contribution:
- Systematic comparison of classical FEM+AD vs PINNs
- Modern implementation using automatic differentiation frameworks
- Demonstration that classical methods often outperform trendy ML approaches

### 🚫 **Why This Won't Scale to Industry Problems**

#### Technical Limitations:
1. **Memory explosion** for large parameter spaces
2. **Gradient computation cost** scales linearly with parameters
3. **No uncertainty quantification**
4. **Requires smooth, well-conditioned problems**

#### Real-world Constraints:
1. **Noisy measurements** 
2. **Model uncertainty**
3. **Real-time decision making**
4. **Distributed/parallel computation needs**

#### Example: Geosteering Applications
- **Parameters**: 10⁶ (permeability/porosity fields)
- **Elements**: 10⁸ (reservoir scale)
- **Time constraint**: Minutes while drilling
- **Uncertainty**: Critical for decision making

**DFEM approach would be computationally prohibitive.**

## Recommendations

### For the Authors:
1. **Honest positioning** - This is "FEM+AD benchmarking against PINNs", not a novel method
2. **Scale analysis** - Demonstrate computational complexity and memory scaling
3. **Industry-scale examples** - Test on realistic problem sizes
4. **Uncertainty quantification** - Essential for practical inverse problems
5. **Comparison with established methods** - EnKF, adjoint methods, etc.

### For the Community:
1. **Value classical methods** - Don't assume ML approaches are always better
2. **Scale matters** - Toy examples don't predict real-world performance
3. **Consider total workflow** - Inverse problems need uncertainty quantification
4. **Historical awareness** - Many "novel" ideas have been explored before

## Conclusion

This paper demonstrates that **classical FEM with automatic differentiation can outperform PINNs** on small-scale inverse problems. However:

- **Novelty claims are overstated** - FEM+AD is decades old
- **Examples are too small** to validate scalability claims  
- **Industry applicability is questionable** due to computational scaling
- **Missing comparisons** with established inverse problem methods

### Bottom Line:
**Useful benchmarking study showing classical methods > trendy ML, but not the breakthrough claimed.**

The real lesson: **Before chasing the latest ML trends, check if classical methods with modern implementation tools can solve your problem better.**

---

## References for Historical Context

- Dixon, L.C.W. et al. (1988). "Finite element optimization: The use of structured automatic differentiation"
- Griewank, A. & Walther, A. (2008). "Evaluating derivatives: principles and techniques of algorithmic differentiation"
- Bartholomew-Biggs, M. et al. (2000). "Automatic differentiation of algorithms"
- MOOSE framework documentation on Automatic Differentiation (2018+)
- Various papers on adjoint methods and EnKF for inverse problems (1990s-2000s)

## how it relates to LSFEM
- For Maxwell equations, LSFEM is not good. Needs to find a variational form of the physical problem. Please refer Jin's book, Chapter 6.

*This review aims to provide balanced perspective on the paper's contributions while highlighting limitations and historical context often missing from the ML-focused literature.*
