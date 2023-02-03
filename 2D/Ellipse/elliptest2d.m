% This script is a test for FEM solver of 2d stable elliptic equation:
%       -div(c*grad(u))+d·grad(u)+k*u=f on \Omega, 
% ( c>0 is a 2x2 positive-definite matrix or a positive scalar, d is a 1x2 vector )
% with three boundary conditions: 
% 1. Dirichlet: u = gD
% 2. Neumann: c*grad(u)·n = gN
% 3. Robin: c*grad(u)·n + gR1*u = gR2
% Remark: the mesh data is got from iFEM software package

%% This script only solves elliptic equation but doesn't computes error

close all;
clear;
% clc;

tic;

elemtype = 3; % 3 for triangle mesh and 4 for quadrilateral mesh
femtype.type = 'pgm'; %PG method
%femtype.type = 'ipdg'; %IPDG method

trbasis_type = 'l'; trbasis_order = 2;
tsbasis_type = 'l'; tsbasis_order = 2;
Gauss_order = 7;
errtype = 'l2';
plottype = '2'; % '1':draw mesh; '2':plot solution;
                % '3':stiffness matrix; '4':cpu time
plotdx_order = 1;
plotdy_order = 0;
upperbd = 1e2;  
refine = 4;

%% example: a polynomial case, \Omega=(-1,1)x(-1,1)
switch elemtype
    case 3
        %[node,elem] = squaremesh([-1,1,-1,1],2); % structure mesh
        [node,elem] = regpolygon(6,1);  % polygon mesh
    case 4
        [node,elem] = squarequadmesh([-1,1,-1,1],2);
end
mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
bdFlag = setboundary(node,elem,'Dirichlet');
bdinfo = genMeshbdinfo(elem,bdFlag);
coeff.c = 10;
coeff.d = [1,0];
coeff.k = -10;
f = @(x,y) 1;
bdinfo.function = cell(1,4);
bdinfo.function{1} = @(x,y) 0+x-x;

for i = 1:refine
    % generate mesh and bdinfo
    switch elemtype
        case 3
            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        case 4
            [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
    end
    mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
    bdinfo = genMeshbdinfo(elem,bdFlag,bdinfo);
    %      if elemtype == 4   % counterclockwise
    %          t = mesh.T(2,:); mesh.T(2,:) = mesh.T(4,:); mesh.T(4,:) = t;
    %      end
end

%% FEM main routine
[uh,~,basis,A,F] = lineEllip2(mesh,coeff,f,bdinfo,trbasis_type,...
    trbasis_order,tsbasis_type,tsbasis_order,Gauss_order,femtype);

fprintf('dof = %d \n',basis.trial.Nb);

toc;

if contains(plottype,'1') % draw mesh
    figure(1);
    showmesh(node,elem);
    %findnode(basis.test.Pb');
    %findelem(node,elem);
end

if contains(plottype,'2') % plot fem solution
    s = num2str(plotdx_order); t = num2str(plotdy_order);
    st = num2str(plotdx_order + plotdy_order);
    str1 = strcat(['FEM linear interpolation solution $I_h(d^',st,'u_h',...
                    '/dx^',s,'dy^',t,')$']);
    str2 = strcat(['Exact solution $I_h(d^',st,'u/dx^',s,'dy^',t,')$']);
    figure(2);
    set(gcf,'unit','normalized','position',[0.2, 0.2, 0.5, 0.7]);
    if basis.trial.type == 'l' && basis.test.type == 'l' && ...
       basis.trial.order == 1 && basis.test.order == 1
        Iuh = uh;
    else
        [Iuh, ~] = utriFEM2(mesh,basis,uh,mesh.P(1,:)',mesh.P(2,:)',...
            plotdx_order,plotdy_order,femtype,upperbd);
    end
    p0 = mesh.P;
    t0 = [mesh.T;ones(1,length(mesh.T(1,:)))];
    if isreal(Iuh)
        pdesurf(p0,t0,Iuh);
    else
        pdesurf(p0,t0,abs(Iuh));
        %pdesurf(p0,t0,-angle(Iuh));
    end
    view(0,90);
    colorbar;
    colormap parula;
    title(str1,'Interpreter','latex');
    set(gca,'fontsize',18);
    grid on;

end

if contains(plottype,'3') % plot stiffness matrix
    figure(3)
    spy(A);
    title("Nonzero element in stiffness matrix");
    set(gca,'fontsize',18);
end

if contains(plottype,'4') % plot cpu time during each refinement
    figure(4)
    plot(1:1:length(T),T,'-o','linewidth',1);
    title("CPU time");
    set(gca,'fontsize',18);
end

 