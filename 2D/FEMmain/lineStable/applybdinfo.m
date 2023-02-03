function [A,F] = applybdinfo(A,F,mesh,trPT,tsPT,fembdinfo)
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
    gD = fembdinfo.function{1};
    bdidx = zeros(Ndofts,1);
    bdidx(Dnodetr) = 1;
    Tbd = spdiags(bdidx,0,Ndofts,Ndoftr);
    T1 = spdiags(1-bdidx,0,Ndofts,Ndofts);
%    T2 = spdiags(1-bdidx,0,Ndoftr,Ndoftr);
    A = T1*A + Tbd;
%    A0 = A;
%    A = T1*A*T2 + Tbd;
%    dA = A0 - A;
end
clear A0;

%% modify load vector
u = zeros(Ndoftr,1);
% Dirichlet boundary condition
%if ~isPureNeumann && ~isempty(fixedNode) && ~isempty(pde.g_D)
if ~isempty(Dnodets)
    pxy = tsPT.Pb(:,Dnodets);
    u(Dnodetr) = gD(pxy(1,:)',pxy(2,:)');
%    F = F - dA*u;
end
if ~isempty(Dnodets) % non-empty Dirichlet boundary condition
    F(Dnodets) = u(Dnodetr);
end


end

