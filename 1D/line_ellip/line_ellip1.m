function [uh,ufem,basis,A,F] = line_ellip1(mesh,coeff,f,bdinfo,...
                           trbasis_type,trbasis_order,tsbasis_type,tsbasis_order,...
                           Gauss_nodes,femtype)
% A FEM solver for 1d stable elliptic equation -(cu')'+du'+ku=f on (a,b),c>0,
% with three boundary conditions: 
% 1. Dirichlet: u=g
% 2. Neumann: u'=r
% 3. Robin: u'+pu=q
% it returns a FEM function uh as the solution
%
% Parameter specification:
% mesh: a struct consisting of mesh info, default: uniform mesh     
% mesh.type: 1 for uniform mesh, 2 for self adaptive mesh, 
%            3 for user-input mesh 
% mesh.left: the left node of the interval, i.e., a 
% mesh.right: the right node of the interval, i.e. b
% mesh.P & mesh.T: the jth column of the matrix P stores the coordinates of  
%                  the jth mesh node; the nth column of the matrix T stores 
%                  the global node indices of the mesh nodes of the nth mesh element
% mesh.N: number of mesh elements
% mesh.Nm: number of mesh nodes
%
% basis: a struct consisting of basis function info, both trial and test
%        function, default: Lagrange linear basis
% basis.trial/test.type: 'l' for Lagrange basis funciton, 
%                        'h' for Hermite basis function
% basis.trial/test.order: a positive integer representing the order
%                         of basis function, less than 7
% basis.trial/test.reffun: a 6x2 cell consisting of all the local basis functions on
%                          reference element, the ith row stores the i-1
%                          derivative of the basis function, the 1st column
%                          is the left nodal basis function and the 2nd the
%                          right one
% basis.trial/test.Pb/Tb: Pb and Tb matrix similar to P and T but storing 
%                         FEM nodes infomation
% reffun编号规则：先从左到右排函数值，再排一阶导数值，再排二阶导数值，以此类推
% basis.trial/test.Nb: number of global basis functions, i.e. number of unknowns
% basis.trial/test.Nlb: number of local basis functions
% 
% bdinfo: shows the boundary condition
% bdinfo.nodes: a 2x2 matrix, the 1st row stores the type of
%               the boundary points (1: Dirichlet, 2: Neumann, 3: Robin), 
%               the 2nd row stores the global index of the boundary points
% bdinfo.value: a 1/2x2 matrix that stores the coefficient value of 
%               boundary condition. For Dirichlet b.c., the corresponding  
%               column is g(a)/g(b); for Neumann b.c., it's r(a)/r(b);
%               for Robin b.c., the 1st row stores p(a)/p(b), the 2nd row 
%               stores q(a)/q(b)
% femtype.type: 'pgm': Petrov-Galerkin method; 
%               'ipdg': IPDG method;
% femtype.ipdg.gamma0, femtype.ipdg.gamma1, femtype.ipdg.alpha, femtype.ipdg.beta:
%               if femtype.type == 'ipdg', the user needs to input the
%               above extra parameters which are vectors with the length
%               the same as mesh.P
% femtype.ipdg.theta: 1 for SIPG; -1 for NIPG; 0 for IIPG

% Example:  
%  coeff.c = @(x) exp(-x);
%  coeff.d = @(x) 0;
%  coeff.k = @(x) 1;
%  f = @(x) (exp(x)*sin(x)+sin(x)-cos(x));
%  bdinfo.flag = [1, 1];
%  bdinfo.value = [0, exp(pi/2)];
%  [uh,mesh,basis] = stab_ellip1(0,1,coeff,f,bdinfo,1,8,'l',1,'l',1,4);

basis.trial.type = trbasis_type;
basis.trial.order = trbasis_order;
basis.test.type = tsbasis_type;
basis.test.order = tsbasis_order;
basis = genbasis1(mesh,basis,femtype); 
[A,F] = assemb_line_ellip1(mesh,basis,coeff,f,bdinfo,Gauss_nodes,femtype);
uh = A\F;  % uh = linear_solver(A,b);
%uh = gmres(A,F,2000);
ufem = @(x) gen_ufem1(mesh,basis,uh,x,femtype);  % generate FEM solution;
end

