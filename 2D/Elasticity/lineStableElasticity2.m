function [uh,ufem,basis,A,F] = lineStableElasticity2(mesh,coeff,f1,f2,bdinfo,trbasis_type,trbasis_order,...
                                tsbasis_type,tsbasis_order,Gauss_order,femtype)
% A FEM solver for 2d linear stable elasticity equations:
%               -div(sigma(u)) = f on \Omega, 
% ( sigma is related to parameter lambda and mu )
% with three boundary conditions: 
% 1. Dirichlet: u = gD
% 2. Neumann: c*grad(u)·n = gN
% 3. Robin: c*grad(u)·n + gR1*u = gR2
% it returns a FEM function uh as the solution
%
% Parameter specification:
% mesh: a struct consisting of mesh info, default: uniform mesh     
% mesh.type: 1.1 for user-input triangle mesh, 1.2 for user-input quadrilateral mesh,
%            2 for self adaptive mesh,  
% mesh.P & mesh.T: the jth column of the matrix P stores the coordinates of  
%                  the jth mesh node; the nth column of the matrix T stores 
%                  the global node indices of the mesh nodes of the nth mesh element
% mesh.N: number of mesh elements
% mesh.Nm: number of mesh nodes
%
% basis: a struct consisting of basis function info, both trial and test
%        function, default: Lagrange linear basis
% basis.trial/test.type: 'l' for Lagrange basis funciton, 
%                        'h' for Hermite basis function,
%                        'cr' for Crouzeix-Raviart element,
%                        'a' for Argyris element,
%                        'b' for Bell element,
%                        'hct' for Hsich-Clough-Tocher element,
%                        'rt' for Raviart-Thomas element
% basis.trial/test.order: a positive integer representing the order
%                         of basis function, less than 8
% basis.trial/test.reffun: a 7 x Nlb x (num_of_derivatives+1) cell consisting of 
%                        all the local basis functions on reference element
%                        the ith row stores the functions on reference element,
%                        the column is corresponding to the nodal basis  
% basis.trial/test.Pb/Tb: Pb and Tb matrix similar to P and T but storing 
%                         FEM nodes infomation
% reffun编号规则：先排函数值，再排一阶导数值，再排二阶导数值，以此类推
% basis.trial/test.Nb: number of global basis functions, i.e. number of unknowns
% basis.trial/test.Nlb: number of local basis functions
% 
% bdinfo: shows the boundary condition of mesh data
% bdinfo.nodes: a (2 x num_of_boundarynodes) matrix 
%               bdinfo.nodes(1,k) is the type of the kth boundary mesh
%                                 node (1: Dirichlet, 2: Neumann, 3: Robin);
%               bdinfo.nodes(2,k) is the global node index of the kth
%                                 boundary mesh node.
% bdinfo.edges: a (4 x num_of_boundaryedges) matrix
%               bdinfo.edges(1,k) is the type of the kth boundary edge ek:
%                                 Dirichlet (1), Neumann (2), Robin (3);
%               bdinfo.edges(2,k) is the index of the element which
%                                 contains the kth boundary edge ek.
%               bdinfo.edges(3,k) is the global node index of the 1st end
%                                 node of the kth boundary edge ek.
%               bdinfo.edges(4,k) is the global node index of the second
%                                 end node of the kth boundary edge ek.
% bdinfo.function: a 1x4 cell that retores function gD, gN, gR1, gR2 for 
%                  Dirichlet, Neumann and Robin condition respectively
%               
% femtype.type: 'pgm': Petrov-Galerkin method; 
%               'ipdg': IPDG method;
% femtype.ipdg.gamma0, femtype.ipdg.gamma1, femtype.ipdg.alpha, femtype.ipdg.beta:
%               if femtype.type == 'ipdg', the user needs to input the
%               above extra parameters which are vectors with the length
%               the same as mesh.P
% femtype.ipdg.theta: 1 for SIPG; -1 for NIPG; 0 for IIPG

elemtype = size(mesh.T,1);
%% generate trial and test basis
switch elemtype
    case 3  % triangular mesh
        basis.trial = genTribasis2(trbasis_type,trbasis_order); 
        basis.test = genTribasis2(tsbasis_type,tsbasis_order);
    case 4  % quadrilateral mesh
        basis.trial = genQuadbasis2(trbasis_type,trbasis_order);
        basis.test = genQuadbasis2(tsbasis_type,tsbasis_order);
end

%% generate Pb and Tb information matrix
trPT = genPbTb(mesh,trbasis_type,trbasis_order,femtype);
tsPT = genPbTb(mesh,tsbasis_type,tsbasis_order,femtype);
basis.trial.Pb = trPT.Pb; basis.trial.Tb = trPT.Tb;
basis.test.Pb = tsPT.Pb; basis.test.Tb = tsPT.Tb;
basis.trial.Nb = size(trPT.Pb,2); basis.trial.Nlb = size(trPT.Tb,1);
basis.test.Nb = size(tsPT.Pb,2); basis.test.Nlb = size(tsPT.Tb,1);

%% assemble stiffness matrix and load vector
switch elemtype
    %% triangular mesh
    case 3  
        A1 = assembTristiff(mesh,basis,coeff.l,1,0,1,0,Gauss_order,femtype);
        A2 = assembTristiff(mesh,basis,coeff.m,1,0,1,0,Gauss_order,femtype);
        A3 = assembTristiff(mesh,basis,coeff.m,0,1,0,1,Gauss_order,femtype);
        A4 = assembTristiff(mesh,basis,coeff.l,1,0,0,1,Gauss_order,femtype);
        A5 = assembTristiff(mesh,basis,coeff.m,0,1,1,0,Gauss_order,femtype);
        A6 = assembTristiff(mesh,basis,coeff.l,0,1,1,0,Gauss_order,femtype);
        A7 = assembTristiff(mesh,basis,coeff.m,1,0,0,1,Gauss_order,femtype);
        A8 = assembTristiff(mesh,basis,coeff.l,0,1,0,1,Gauss_order,femtype);
        A = [A1 + 2*A2 + A3, A4 + A5; A6 + A7, A8 + 2*A3 + A2];
        clear A1 A2 A3 A4 A5 A6 A7 A8;
        F1 = assembTriload(mesh,basis,f1,0,0,Gauss_order,femtype);
        F2 = assembTriload(mesh,basis,f2,0,0,Gauss_order,femtype);
        F = [F1; F2];
        clear F1 F2;
    %% quadrilateral mesh
    case 4  
        
end


%% add boundary condition
fembdinfo.trial = genFEMbdinfo(bdinfo,trbasis_type,trbasis_order,trPT,elemtype,femtype);
fembdinfo.test = genFEMbdinfo(bdinfo,tsbasis_type,tsbasis_order,tsPT,elemtype,femtype);
fembdinfo.function = bdinfo.function;
[A,F] = applyEqsbdinfo(A,F,mesh,trPT,tsPT,fembdinfo);

%% solve pde and get FEM solution
uh = A\F;  % uh = FEMsolver(A,b);
ufem = @utriFEM2;  % generate FEM solution;

end

% subroutine to apply boundary condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,F] = applyEqsbdinfo(A,F,mesh,trPT,tsPT,fembdinfo)
% apply boundary information to modify stiffness matrix and load vector for 2d 

% fembdinfo: a struct consists the following member
% fembdinfo.trial/test.nodes: a (2 x num_of_boundarynodes) matrix 
%               nodes(1,k) is the type of the kth boundary finite
%                          element node (1: Dirichlet, 2: Neumann, 3: Robin);
%               nodes(2,k) is the global node index of the kth
%                          boundary finite element node.
% fembdinfo.trial/test.edges: a (4 x num_of_boundaryedges) matrix
%               edges(1,k) is the type of the kth boundary edge ek:
%                          Dirichlet (1), Neumann (2), Robin (3);
%               edges(2,k) is the index of the element which
%                          contains the kth boundary edge ek.
%               edges(3,k) is the global node index of the 1st end
%                          node of the kth boundary edge ek.
%               edges(4,k) is the global node index of the second
%                          end node of the kth boundary edge ek.
% fembdinfo.function: a 1x4 cell that retores function gD, gN, gR1, gR2
%               for Dirichlet, Neumann and Robin condition respectively


[Ndofts, Ndoftr] = size(A);

% Dirichlet boundary node
Dnodetr = fembdinfo.trial.nodes(2,fembdinfo.trial.nodes(1,:)==1);
Dnodets = fembdinfo.test.nodes(2,fembdinfo.test.nodes(1,:)==1);
fixedNode = Dnodetr; 
% Neumann boundary node
Nnode = fembdinfo.trial.nodes(2,fembdinfo.trial.nodes(1,:)==2);
% Robin boundary node
Rnode = fembdinfo.trial.nodes(2,fembdinfo.trial.nodes(1,:)==3);

freeNode = 1:Ndoftr;
freeNode(fixedNode)=[];

%% modify stiffness matrix
% Dirichlet boundary condition
% A(Dnode,Dnode)=I, A(Dnode,~Dnode)=0, A(~Dnode,Dnode)=0
if ~isempty(Dnodetr)
    g1D = fembdinfo.function{1,1};
    g2D = fembdinfo.function{2,1};
    bdidx = zeros(Ndofts/2,1);
    bdidx(Dnodetr) = 1;
    Tbd = spdiags([bdidx;bdidx],0,Ndofts,Ndoftr);
    T1 = spdiags([1-bdidx;1-bdidx],0,Ndofts,Ndofts);
    A = T1*A + Tbd;
end
clear A0;

%% modify load vector
u = zeros(Ndoftr,1);
% Dirichlet boundary condition
%if ~isPureNeumann && ~isempty(fixedNode) && ~isempty(pde.g_D)
if ~isempty(Dnodets)
    pxy = tsPT.Pb(:,Dnodets);
    u(Dnodetr) = g1D(pxy(1,:)',pxy(2,:)');
    u(Dnodetr + Ndoftr/2) = g2D(pxy(1,:)',pxy(2,:)');
%    F = F - dA*u;
end
if ~isempty(Dnodets) % non-empty Dirichlet boundary condition
    F(Dnodets) = u(Dnodetr);
    F(Dnodets + Ndofts/2) = u(Dnodetr + Ndoftr/2);
end


end


