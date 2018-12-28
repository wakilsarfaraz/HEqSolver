# HEqSolver
HEqSolver is a Finite Element Solver written in R for solving the diffusion equation. HEqSolver employes piecewise linear basis functions on uniform triangulation. It solves the diffusion equation on a two dimensional rectangular domain with homogeneous boundary conditions of Dirichlet type. We illustrate HEqSolver through an example of the diffusion equation posed on a rectangular domain in which we present both the relevant theory of the finite element method a long with the time-stepping scheme. All this is illustrated through the presentation of the mathematical theory of the finite element method on an example, for which the solution is obtained by HEqSolver. 
## Finite element method overview
Consider a two dimensional open bounded [convex](https://en.wikipedia.org/wiki/Convex_set) domain denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega\subset\mathbb{R}^2"/> on the positive real <img src="https://latex.codecogs.com/svg.latex?\Large&space;(x,y)\in\mathbb{R}^2"/> plane. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,t)"/> satisfy the probelm <img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{\partial\,u}{\partial\,t}=D\Delta\,u(x,y)"/> for <img src="https://latex.codecogs.com/svg.latex?\Large&space;\!\,x,y\in\Omega,\,\,\,t>0"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,t)=0"/> for <img src="https://latex.codecogs.com/svg.latex?\Large&space;\quad\,x,y\in\partial\Omega"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;t>0"/>. Initial conditions are prescribed such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,0)=u_0(x,y)"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;u_0(x,y)\in\mathbb{L}_2"/>. For solving time-dependent problems through numerical methods two types of approaches are possible. The two approaches mainly depend on whether one solves in time first then in space or vice versa. We adopt the approach in which we apply Galerkin finite element method to solve the equation in space first and then we proceed with time-stepping using the implicit Euler scheme. Let the solution and test spaces be denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}=\{w\in\,H^1:\|w(x,y)\|^2<\infty,\,\|\nabla\,w(x,y)\|^2<\infty\}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}=\{w\in\,H^1:\,w(x,y)=0,\;\,x,y\in\partial\Omega\}"/> respectively, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are a special class of [Hilbert spaces](https://en.wikipedia.org/wiki/Hilbert_space) known as [Sobelov spaces](https://en.wikipedia.org/wiki/Hilbert_space). Assuming that the problem satisfies all the conditions stated by [Lax-Milgram Theorem](http://mathworld.wolfram.com/Lax-MilgramTheorem.html), we multiply both sides of the equation by a test function <img src="https://latex.codecogs.com/svg.latex?\Large&space;w(x,y)\in\,H^1_{\Omega_0}"/> and integrate by parts with application of [Green's Theorem](https://en.wikipedia.org/wiki/Green%27s_theorem), then the weak formulation of the problem is to find <img src="https://latex.codecogs.com/svg.latex?\Large&space;u\in\,H^1_{\Omega}"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega}\nabla\!u\cdot\nabla\!w\,d\Omega=\int_{\Omega}fw\,d\Omega"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;w\in\,H^1_{\Omega_0}"/>. Given that the problem is [well-posed](https://en.wikipedia.org/wiki/Well-posed_problem) in Hilbert space of functions of which both <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are subsets, then we can define finite dimensional solution and test spaces namely <img src="https://latex.codecogs.com/svg.latex?\Large&space;V^h\subset\,H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;V_0^h\subset\,H^1_{\Omega_0}"/> respectively. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> be a triangulated approximation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, containing <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> discrete number of nodes and <img src="https://latex.codecogs.com/svg.latex?\Large&space;K"/> triangles, then using the finite dimentional solution and test spaces, the finite element formulation of the problem is to find <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h\in\,V^h"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega^h}\nabla\!u^h\cdot\nabla\!w^h\,d\Omega=\int_{\Omega^h}fw^h\,d\Omega"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;\,w^h\in\,V_0^h"/>.  We expand <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h"/> in terms of its finite i.e. <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> dimensional [basis functions](https://en.wikipedia.org/wiki/Basis_function) in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h=\sum_{i=1}^N\,U_i\phi_i"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h=\sum_{j=1}^N\,\phi_j"/>. Substituting these in the finite element formulation we obtain a discrete set of linear equations of the from <img src="https://latex.codecogs.com/svg.latex?\Large&space;S\,U=L"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> is the stiffness matrix, <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> is the vector of unknowns containing the approximate solution values at each node and <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> is the right hand-side load vector. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> are given by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\int_{\Omega^h}\,\nabla\phi_i\cdot\nabla\phi_j\,d\Omega^h"/>, <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[F\]_{j}=\int_{\Omega^h}\,f\phi_j\,d\Omega^h"/>. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> are the finite element approximate solution values to be found from the linear system. 
For code implementation we start with a scheme to descritise  <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>. This is achieved by creating an <img src="https://latex.codecogs.com/svg.latex?\Large&space;N\times\,N"/> grid points within <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, which provides a uniform quadrilateral mesh. Then we triangulate the quadrilateral mesh by drawing a line of slope -1 diangonally through each square in order to divide every quadrilateral element into two adjacent triangles, which leads to the construction of a uniform triangulated domain <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/>.  Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathcal{T}"/> denote the uniform triangulation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, consisting of triangles. The algorithm is programmed to compute locally on each triangle <img src="https://latex.codecogs.com/svg.latex?\Large&space;K\in\mathcal{T}"/>, the entries of the matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and the load vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/>. The entries of the local <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times3"/> matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and  <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times1"/> vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> are constructed by the formulea <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\sum_{K}\int_{\Omega^h}\nabla\phi_i\cdot\nabla\phi_j\,d\Omega"/>  and <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[L\]_j=\sum_K\int_{\Omega^h}f\phi_j\,d\Omega"/>.
