function [uh,ufem,basis,A,F] = lineEllip2(mesh,coeff,f,bdinfo,trbasis_type,trbasis_order,...
                                tsbasis_type,tsbasis_order,Gauss_order,femtype)
% A FEM solver for 2d stable elliptic equation:
%               -div(c*grad(u)) + d·grad(u) + k*u = f on \Omega, 
% ( c>0 is a 2x2 matrix, d is a 1x2 vector )
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

%% assemble stiffness matrix and load vector, trial basis = test basis
switch elemtype
    %% triangular mesh
    case 3  
        % Part 1: assemble -div(c*grad(u))
        if iscell(coeff.c)  % c is a 2x2 matrix with each element a function
            if size(coeff.c,2) == 1  % c is a mx1 vector
                coeff.c = coeff.c';
            end
            switch size(coeff.c,1)
                case 1  % c is a 1xm vector, m=1,2,3,4
                    switch size(coeff.c,2)
                        case 1  % c is an identity plus a scalar function
                            A11 = assembTristiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c{1},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 2  % c is diagonal matrix function
                            A11 = assembTristiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c{2},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 3  % c is symmetric matrix function
                            A11 = assembTristiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A12 = assembTristiff(mesh,basis,coeff.c{3},0,1,1,0,Gauss_order,femtype);
                            A13 = assembTristiff(mesh,basis,coeff.c{3},1,0,0,1,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c{2},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                        case 4  % c is non-symmetric matrix function
                            A11 = assembTristiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A12 = assembTristiff(mesh,basis,coeff.c{3},0,1,1,0,Gauss_order,femtype);
                            A13 = assembTristiff(mesh,basis,coeff.c{4},1,0,0,1,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c{2},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                    end

                case 2  % c is a 2x2 matrix
                    A11 = assembTristiff(mesh,basis,coeff.c{1,1},1,0,1,0,Gauss_order,femtype);
                    A12 = assembTristiff(mesh,basis,coeff.c{1,2},0,1,1,0,Gauss_order,femtype);
                    A13 = assembTristiff(mesh,basis,coeff.c{2,1},1,0,0,1,Gauss_order,femtype);
                    A14 = assembTristiff(mesh,basis,coeff.c{2,2},0,1,0,1,Gauss_order,femtype);
                    A1 = A11 + A12 + A13 + A14;
                    clear A11 A12 A13 A14;
            end
        elseif isnumeric(coeff.c) % c is a vector or matrix with each element a scalar
            if size(coeff.c,2) == 1  % c is a mx1 vector
                coeff.c = coeff.c';
            end
            switch size(coeff.c,1)
                case 1  % c is a 1xm vector, m=1,2,3,4
                    switch size(coeff.c,2)
                        case 1  % c is an identity plus a scalar
                            A11 = assembTristiff(mesh,basis,coeff.c,1,0,1,0,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c,0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 2  % c is diagonal matrix
                            A11 = assembTristiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c(2),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 3  % c is symmetric matrix
                            A11 = assembTristiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A12 = assembTristiff(mesh,basis,coeff.c(3),0,1,1,0,Gauss_order,femtype);
                            A13 = assembTristiff(mesh,basis,coeff.c(3),1,0,0,1,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c(2),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                        case 4  % c is non-symmetric matrix
                            A11 = assembTristiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A12 = assembTristiff(mesh,basis,coeff.c(3),0,1,1,0,Gauss_order,femtype);
                            A13 = assembTristiff(mesh,basis,coeff.c(4),1,0,0,1,Gauss_order,femtype);
                            A14 = assembTristiff(mesh,basis,coeff.c(2),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                    end

                case 2  % c is a 2x2 matrix
                    A11 = assembTristiff(mesh,basis,coeff.c(1,1),1,0,1,0,Gauss_order,femtype);
                    A12 = assembTristiff(mesh,basis,coeff.c(1,2),0,1,1,0,Gauss_order,femtype);
                    A13 = assembTristiff(mesh,basis,coeff.c(2,1),1,0,0,1,Gauss_order,femtype);
                    A14 = assembTristiff(mesh,basis,coeff.c(2,2),0,1,0,1,Gauss_order,femtype);
                    A1 = A11 + A12 + A13 + A14;
                    clear A11 A12 A13 A14;
            end
        else  % c is a scalar function
            A11 = assembTristiff(mesh,basis,coeff.c,1,0,1,0,Gauss_order,femtype);
            A12 = assembTristiff(mesh,basis,coeff.c,0,1,0,1,Gauss_order,femtype);
            A1 = A11 + A12;
            clear A11 A12;
        end

        % Part 2: assemble d·grad(u)
        if iscell(coeff.d) % d is a 1x2 vector with each element a function
            A21 = assembTristiff(mesh,basis,coeff.d{1},1,0,0,0,Gauss_order,femtype);
            A22 = assembTristiff(mesh,basis,coeff.d{2},0,1,0,0,Gauss_order,femtype);
            A2 = A21 + A22;
            clear A21 A22;
        else % d is a 1x2 vector with each element a constant
            A21 = assembTristiff(mesh,basis,coeff.d(1),1,0,0,0,Gauss_order,femtype);
            A22 = assembTristiff(mesh,basis,coeff.d(2),0,1,0,0,Gauss_order,femtype);
            A2 = A21 + A22;
            clear A21 A22;
        end

        % Part 3: assemble k*u
        A3 = assembTristiff(mesh,basis,coeff.k,0,0,0,0,Gauss_order,femtype);
        % add three part on
        A = A1 + A2 + A3;
        clear A1 A2 A3

        % Part 4: assemble f
        F = assembTriload(mesh,basis,f,0,0,Gauss_order,femtype);

        
    %% quadrilateral mesh
    case 4  
        % Part 1: assemble -div(c*grad(u))
        if iscell(coeff.c)  % c is a 2x2 matrix with each element a function
            if size(coeff.c,2) == 1  % c is a mx1 vector
                coeff.c = coeff.c';
            end
            switch size(coeff.c,1)
                case 1  % c is a 1xm vector, m=1,2,3,4
                    switch size(coeff.c,2)
                        case 1  % c is an identity plus a scalar function
                            A11 = assembQuadstiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c{1},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 2  % c is diagonal matrix function
                            A11 = assembQuadstiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c{2},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 3  % c is symmetric matrix function
                            A11 = assembQuadstiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A12 = assembQuadstiff(mesh,basis,coeff.c{3},0,1,1,0,Gauss_order,femtype);
                            A13 = assembQuadstiff(mesh,basis,coeff.c{3},1,0,0,1,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c{2},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                        case 4  % c is non-symmetric matrix function
                            A11 = assembQuadstiff(mesh,basis,coeff.c{1},1,0,1,0,Gauss_order,femtype);
                            A12 = assembQuadstiff(mesh,basis,coeff.c{3},0,1,1,0,Gauss_order,femtype);
                            A13 = assembQuadstiff(mesh,basis,coeff.c{4},1,0,0,1,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c{2},0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                    end

                case 2  % c is a 2x2 matrix
                    A11 = assembQuadstiff(mesh,basis,coeff.c{1,1},1,0,1,0,Gauss_order,femtype);
                    A12 = assembQuadstiff(mesh,basis,coeff.c{1,2},0,1,1,0,Gauss_order,femtype);
                    A13 = assembQuadstiff(mesh,basis,coeff.c{2,1},1,0,0,1,Gauss_order,femtype);
                    A14 = assembQuadstiff(mesh,basis,coeff.c{2,2},0,1,0,1,Gauss_order,femtype);
                    A1 = A11 + A12 + A13 + A14;
                    clear A11 A12 A13 A14;
            end
        elseif isnumeric(coeff.c) % c is a vector or matrix with each element a scalar
            if size(coeff.c,2) == 1  % c is a mx1 vector
                coeff.c = coeff.c';
            end
            switch size(coeff.c,1)
                case 1  % c is a 1xm vector, m=1,2,3,4
                    switch size(coeff.c,2)
                        case 1  % c is an identity plus a scalar
                            A11 = assembQuadstiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c(1),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 2  % c is diagonal matrix
                            A11 = assembQuadstiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c(2),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A14;
                            clear A11 A14;
                        case 3  % c is symmetric matrix
                            A11 = assembQuadstiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A12 = assembQuadstiff(mesh,basis,coeff.c(3),0,1,1,0,Gauss_order,femtype);
                            A13 = assembQuadstiff(mesh,basis,coeff.c(3),1,0,0,1,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c(2),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                        case 4  % c is non-symmetric matrix
                            A11 = assembQuadstiff(mesh,basis,coeff.c(1),1,0,1,0,Gauss_order,femtype);
                            A12 = assembQuadstiff(mesh,basis,coeff.c(3),0,1,1,0,Gauss_order,femtype);
                            A13 = assembQuadstiff(mesh,basis,coeff.c(4),1,0,0,1,Gauss_order,femtype);
                            A14 = assembQuadstiff(mesh,basis,coeff.c(2),0,1,0,1,Gauss_order,femtype);
                            A1 = A11 + A12 + A13 + A14;
                            clear A11 A12 A13 A14;
                    end

                case 2  % c is a 2x2 matrix
                    A11 = assembQuadstiff(mesh,basis,coeff.c(1,1),1,0,1,0,Gauss_order,femtype);
                    A12 = assembQuadstiff(mesh,basis,coeff.c(1,2),0,1,1,0,Gauss_order,femtype);
                    A13 = assembQuadstiff(mesh,basis,coeff.c(2,1),1,0,0,1,Gauss_order,femtype);
                    A14 = assembQuadstiff(mesh,basis,coeff.c(2,2),0,1,0,1,Gauss_order,femtype);
                    A1 = A11 + A12 + A13 + A14;
                    clear A11 A12 A13 A14;
            end
        else  % c is a scalar function
            A11 = assembQuadstiff(mesh,basis,coeff.c,1,0,1,0,Gauss_order,femtype);
            A12 = assembQuadstiff(mesh,basis,coeff.c,0,1,0,1,Gauss_order,femtype);
            A1 = A11 + A12;
            clear A11 A12;
        end

        % Part 2: assemble d·grad(u)
        if iscell(coeff.d) % d is a 1x2 vector with each element a function
            A21 = assembQuadstiff(mesh,basis,coeff.d{1},1,0,0,0,Gauss_order,femtype);
            A22 = assembQuadstiff(mesh,basis,coeff.d{2},0,1,0,0,Gauss_order,femtype);
            A2 = A21 + A22;
            clear A21 A22;
        else % d is a 1x2 vector with each element a constant
            A21 = assembQuadstiff(mesh,basis,d(1),1,0,0,0,Gauss_order,femtype);
            A22 = assembQuadstiff(mesh,basis,d(2),0,1,0,0,Gauss_order,femtype);
            A2 = A21 + A22;
            clear A21 A22;
        end

        % Part 3: assemble k*u
        A3 = assembQuadstiff(mesh,basis,coeff.k,0,0,0,0,Gauss_order,femtype);
        % add three part on
        A = A1 + A2 + A3;
        clear A1 A2 A3;

        % Part 4: assemble f
        F = assembQuadload(mesh,basis,f,0,0,Gauss_order,femtype);
end


%% add boundary condition
fembdinfo.trial = genFEMbdinfo(bdinfo,trbasis_type,trbasis_order,trPT,elemtype,femtype);
fembdinfo.test = genFEMbdinfo(bdinfo,tsbasis_type,tsbasis_order,tsPT,elemtype,femtype);
fembdinfo.function = bdinfo.function;
[A,F] = applybdinfo(A,F,mesh,trPT,tsPT,fembdinfo);

%% solve pde and get FEM solution
uh = A\F;  % uh = FEMsolver(A,b);
ufem = @utriFEM2;  % generate FEM solution;

end

