# HEqSolver
HEqSolver is a Finite Element Solver written in R for solving the diffusion equation. HEqSolver employes piecewise linear basis functions on uniform triangulation. It solves the diffusion equation on a two dimensional rectangular domain with homogeneous boundary conditions of Dirichlet type. We illustrate HEqSolver through an example of the diffusion equation posed on a rectangular domain in which we present both the relevant theory of the finite element method a long with the time-stepping scheme. All this is illustrated through the presentation of the mathematical theory of the finite element method on an example, for which the solution is obtained by HEqSolver. 
## Finite element method overview
Consider <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega=(0\,\,\,1)\times(0\,\,\,1)\subset\mathbb{R}^2"/>, a rectangular open bounded domain with <img src="https://latex.codecogs.com/svg.latex?\Large&space;\partial\Omega\subset\mathbb{R}^2"/> denoting its boundary. We define the initial conditions on each point in space <img src="https://latex.codecogs.com/svg.latex?\Large&space;(x,y)\in\Omega"/> at time <img src="https://latex.codecogs.com/svg.latex?\Large&space;t=0"/> and hence consider the diffusion equation in two space and one time dimension for <img src="https://latex.codecogs.com/svg.latex?\Large&space;0<t\leq\,T_{max}"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;T_{max}"/> denotes the final time.  Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,t)"/> satisfy the probelm <img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{\partial\,u}{\partial\,t}=D\Delta\,u(x,y)"/> for <img src="https://latex.codecogs.com/svg.latex?\Large&space;\!(x,y)\in\Omega,\,\,\,t\in(0\,\,\,T_{max}\]"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,t)=0"/> for <img src="https://latex.codecogs.com/svg.latex?\Large&space;\quad(x,y)\in\partial\Omega"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;t>0"/>. Here, <img src="https://latex.codecogs.com/svg.latex?\Large&space;D>0"/> denotes the diffusion coefficient. Initial conditions for the problem are prescribed such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,0)=u_0(x,y)"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;u_0(x,y)\in\mathbb{L}_2(\Omega)"/>, which means that it is a square integrable function of <img src="https://latex.codecogs.com/svg.latex?\Large&space;x"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;y"/>. The approximate solution of time-dependent problems through numerical methods can in general be achieved through two different types of approaches. The two approaches mainly depend on whether the scheme is formulated to discretise time first followed by the space discretisation or vice vera. We adopt the approach in which we apply Galerkin finite element method to implement the space discretisation first and then we proceed with time-stepping using the implicit Euler scheme. Consider the diffusion problem to be posed on a space-time cylinder <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega\times(0\,\,\,T_{max}\]"/>, in the form of a sequence of Poisson's equations, that are all well-posed on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, i.e. satisfying all the requirements stated by the [Lax-Milgram Theorem](http://mathworld.wolfram.com/Lax-MilgramTheorem.html) at each time point in the interval <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0\,\,\,T_{max}\]"/>. 
This allows to introduce at each time point <img src="https://latex.codecogs.com/svg.latex?\Large&space;t"/>, the solution and test spaces of functions for the corresponding sequence of Poisson equations to be denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}=\{w\in\,H^1:\|w(x,y)\|^2<\infty,\,\|\nabla\,w(x,y)\|^2<\infty\}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}=\{w\in\,H^1:\,w(x,y)=0,\;\,x,y\in\partial\Omega\}"/> respectively, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are a special class of [Sobelov spaces](https://en.wikipedia.org/wiki/Hilbert_space). Assuming that the diffusion problem is well-posed on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;0<t\leq\,T"/>, then we multiply both sides of the diffusion equation by a test function <img src="https://latex.codecogs.com/svg.latex?\Large&space;w(x,y)\in\,H^1_{\Omega_0}"/> and integrate by parts with application of [Green's Theorem](https://en.wikipedia.org/wiki/Green%27s_theorem). Using the prescribed boundary conditions lead to vanish the boundary terms in the application of integration by parts to the right hand-side of the equation, which results in the weak formulation of the diffusion problem, which is to find some <img src="https://latex.codecogs.com/svg.latex?\Large&space;u\in\,H^1_{\Omega}"/> at each time point <img src="https://latex.codecogs.com/svg.latex?\Large&space;t"/>  such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega}\frac{\partial\,u}{\partial\,t}w\,dxdy=-D\int_{\Omega}\nabla\,u\,\cdot\nabla\,w\,dxdy"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;w\in\,H^1_{\Omega_0}"/>. Given that the problem is [well-posed](https://en.wikipedia.org/wiki/Well-posed_problem) in the continuous space of functions of which both <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are subsets, we define the finite dimensional solution and test spaces namely <img src="https://latex.codecogs.com/svg.latex?\Large&space;V^h\subset\,H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;V_0^h\subset\,H^1_{\Omega_0}"/> respectively for each time point <img src="https://latex.codecogs.com/svg.latex?\Large&space;t"/>. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> be a triangulated approximation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, containing a discrete integer number of nodes and <img src="https://latex.codecogs.com/svg.latex?\Large&space;K"/> triangles, then using the finite dimentional solution and test spaces, the finite element formulation of the problem is to find <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h\in\,V^h"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega^h}\frac{\partial\,u^h}{\partial\,t}w^h\,dxdy=-D\int_{\Omega^h}\nabla\,u^h\cdot\nabla\,w^h\,dxdy"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;\,w^h\in\,V_0^h"/>.  We expand <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h"/> in terms of its finite i.e. <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> dimensional [basis functions](https://en.wikipedia.org/wiki/Basis_function) in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h=\sum_{i=1}^N\,U_i\phi_i"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h=\sum_{j=1}^N\,\phi_j"/>. Substituting these in the finite element formulation we obtain a discrete set of continuous oridnary differential equations in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{d}{dt}\Big(\int_{\Omega}\sum_{i=1}^NU_i(t)\phi_i\sum_{j=1}^N\phi_j\,dxdy\Big)=-D\int_{\Omega}\sum_{i=1}^NU_i(t)\nabla\,u\cdot\sum_{j=1}^N\nabla\,\phi_j\,dxdy"/>. We use the commutative property of summation and integration to write the semi-discrete system of differential equations in matrix notation in the form<img src="https://latex.codecogs.com/svg.latex?\Large&space;M\frac{dU}{dt}=-DSU(t)"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;M"/> respectively denote the stiffness and mass matrices. The vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> contains the unknown approximate solution values at each node in the domain <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> at time <img src="https://latex.codecogs.com/svg.latex?\Large&space;t"/>. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;M"/> are given by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\int_{\Omega^h}\,\nabla\phi_i\cdot\nabla\phi_j\,dxdy"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[M\]_{i,j}=\int_{\Omega^h}\phi_i\phi_j\,dxdy"/> respectively. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> are the finite element approximate solution values to be found from the system of ordinary differential equations. 
In order to solve the semi-discrete problem we introduce a uniform mesh of <img src="https://latex.codecogs.com/svg.latex?\Large&space;T"/> discrete number of time points in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0,\Delta\,t,2\Delta\,t,3\Delta\,t...,T\Delta\,t=T_{max})"/>, which is equally spaced by time step-size <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,t"/>. We use the [Implicit Euler Method](https://en.wikipedia.org/wiki/Backward_Euler_method), to fully discretise the semi-discrete system of equations. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;U^m"/> denote the finite element approximate solution at the time-step <img src="https://latex.codecogs.com/svg.latex?\Large&space;m\Delta\,t"/>, then full discretisation of the problem is achieved by writing <img src="https://latex.codecogs.com/svg.latex?\Large&space;M\frac{U^{m+1}-U^m}{\Delta\,t}=-DSU^{m+1}"/>, which can be rearranged in the form of a system of <img src="https://latex.codecogs.com/svg.latex?\Large&space;T"/> linear equations for each time-step written as <img src="https://latex.codecogs.com/svg.latex?\Large&space;U^{m+1}=\Big(S+D\Delta\,t\,M\Big)^{-1}MU^{m}"/>. 
Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/> be discretised by <img src="https://latex.codecogs.com/svg.latex?\Large&space;2N^2"/> uniform triangles connected through <img src="https://latex.codecogs.com/svg.latex?\Large&space;(N+1)^2"/> nodes.  Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathcal{T}"/> denote the uniform triangulation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, consisting of all triangles. The algorithm is programmed to compute locally on each triangle <img src="https://latex.codecogs.com/svg.latex?\Large&space;K\in\mathcal{T}"/>, the entries of the mass and stiffness matrices namely <img src="https://latex.codecogs.com/svg.latex?\Large&space;M"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/>. The entries of the global mass and stiffness matrices are computed by the local <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times3"/> mass and stiffness matrices namely  <img src="https://latex.codecogs.com/svg.latex?\Large&space;M"/> and  <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/>. The global mass and stiffness matrices are therefore, given by the formulea <img src="https://latex.codecogs.com/svg.latex?\Large&space;M_{i,j}=\sum_{K\in\mathcal{T}}\int_{K}\phi_i\phi_jdxdy"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;S_{i,j}=\sum_{K\in\mathcal{T}}\int_{K}\nabla\phi_i\cdot\nabla\phi_jdxdy"/> respectively. 
## Code documentation
Recall the problem formulated in terms of <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,t)"/> satisfying the the diffusion equation <img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{\partial\,u}{\partial\,t}=\Delta\,u(x,y)"/> on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega=(0\;1)\times(0\;1)"/>. It satisfies homogeneous Dirichlet type boundary conditions for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;t\in(0\,\,\,1\]"/>, which is formally written as <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,t)=0,\;\,(x,y)\in\partial\Omega"/>. The initial conditions are prescribed such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y,0)=1,\;\,(x,y)\in\Omega,\;\,t=0"/>.
We present the code in the form of a serie of enumerated documentations for each block of code:

1. Variable
	``` r
	L = 1
	```
	stores the value for the side length of the two-dimensional square domain <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>.

2. Variables
	``` r
	tm = 1
	dt = 0.001
	```
	store the values respectively for the final time <img src="https://latex.codecogs.com/svg.latex?\Large&space;T_{max}"/> and the 	time step-size <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,t"/>. 

3. The number of time points <img src="https://latex.codecogs.com/svg.latex?\Large&space;T"/> on the time interval <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0\;1\]"/> are stored by
	``` r
	M = tm/dt                                                                                            
	```

4. The number of spatial uniform mesh points that discretises `L` is stored by
	``` r
	N = 50
	```

5. The side length `L` is discretised by `N+1` equally spaced points using   
	``` r
	X = seq(0,L,len=N+1)
	```

6. The time interval is also uniformly discretised by `M+1` points using
	``` r
	T = seq(0,tm,len=M+1)
	```

7. This step is to compute and store the matrices containing the values for the two dimensional spatial coordinates of such 			discretisation, which is achieved by 
	``` r
	m = length(X); n=length(X);
	x = matrix(rep(X,each=n),nrow=n);
	y = matrix(rep(X,m),nrow=n)
	```

8. Reshape the data structure of the spatial coordinates from matrices into vector by 
	``` r
	x = c(x)
	y = c(y)
	```
	where `x` and `y` now store all the values for the x and y coordinates for all the global nodes, therefore, each one is now a 		vector of <img src="https://latex.codecogs.com/svg.latex?\Large&space;(N+1)^2"/> entries. 

9. Store the number of global nodes in the discretised domain by 
	``` r
	GNodes = (N+1)^2
	```

10. Use `U` to store the approximate discrete solution values at each time step and we define the initial state of `U` to be a vector of   	constant values of 1 at each note in the domain and this is achieved by 
	``` r
	U = matrix(1, GNodes, 1)
	```

11. Triangulation of a quadrilateral mesh by drawing a line of slope -1 diangonally through each square leads to the construction of a 		uniform triangulated domain <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/>, which consists of `NumTRI` 	triangles, which is defined in terms of `N` by
	``` r
	NumTRI = 2*N^2
	```

12. The connectivity array `LocNodes` is a matrix that stores the global counting of all the nodes with a local structure in the sense 		that it must have `NumTRI` rows and three columns for each vertix of each triangle. The code to achieve this is
	``` r
	LocNodes = matrix(0,NumTRI,3)
	``` 

13. Triangulation of the domain is obtained such that the vertices of each triangle is locally counted in anti-clockwise orientation and 	such that the local counting of vertices of each triangle fills the global connectivity array `LocNodes` in the correct order. Note 	that we must fill in `4+2=6` vertices for each square, because each square is divided into two triangles by a diagonal line, however, 	each triangle has one distinct but two shared vertices. Therefore, it can be noted that the first three lines inside the nested for loop 	correspond to the three vertices of all the lower triangles on each square and the latter three lines correspond to all the three 	vertices of all the upper triangles. The nested loop that distributes such order of counting within the connectivity array is given 	by 
	``` r
	for (i in 1:N){
		for (j in 1:N){
			    LocNodes[i+2*(j-1)*N,1] = i+(j-1)*(N+1)
	        LocNodes[i+2*(j-1)*N,2] = i+j*(N+1)
	        LocNodes[i+2*(j-1)*N,3] = (i+1)+(j-1)*(N+1)
	        LocNodes[i+N+2*(j-1)*N,1] = i+1+j*(N+1)
	        LocNodes[i+N+2*(j-1)*N,2] = (i+1)+(j-1)*(N+1)
	        LocNodes[i+N+2*(j-1)*N,3] = i+j*(N+1)
	    }
	}
	```

14. Introduce the global stiffness and mass matrices filled with zeros each of size `GNodes x GNodes`  
	``` r
	SP_Stiff <- matrix(0, GNodes, GNodes)
	SP_Mass <- matrix(0, GNodes, GNodes)
	```
15. The next segment of code is a for loop that executes a series of computations on each one of the triangles. It starts with 
	``` r
	for (n in 1:NumTRI)
		{# This is the start of the for loop on all triangles
	```
16. The first part within the loop defines the local position vectors `r1`, `r2`, and `r3` in terms of `x` and `y` coordinate values 		through which we refer to different vertices of each triangle. 
	``` r
	r1 = matrix(c(x[LocNodes[n,1]],y[LocNodes[n,1]]),nrow=2, byrow=FALSE)
	r2 = matrix(c(x[LocNodes[n,2]],y[LocNodes[n,2]]),nrow=2, byrow=FALSE);
	r3 = matrix(c(x[LocNodes[n,3]],y[LocNodes[n,3]]),nrow=2, byrow=FALSE);
	```
17. The second part within the loop defines a <img src="https://latex.codecogs.com/svg.latex?\Large&space;2\times2"/> jacobian matrix `J` 	for the mapping from an arbitrary triangle in <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> to a reference 	triangle that has vertices in order <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0,0),\!(1,0)"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0,1)"/>.
	``` r
	J = matrix(c(r2[1]-r1[1],r2[2]-r1[2]
	    ,r3[1]-r1[1],r3[2]-r1[2]), nrow=2, byrow=TRUE) 
	```
	The Jacobian of the mapping namely `J` serves to reduce the computational cost by a significant amount, particularly due to a 		property of integration for computing the integral of a function on a reference domain with a given mapping between  the arbitrary 	domain and the reference domain. Further details on this topic can be found on [Integral domain transformation](http://www.iue.tuwien.ac.at/phd/nentchev/node58.html). 

	
