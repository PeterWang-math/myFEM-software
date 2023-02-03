% This script is a test for FEM solver of 2d stable elliptic equation:
%       -div(c*grad(u))+d·grad(u)+k*u=f on \Omega, 
% ( c>0 is a 2x2 positive-definite matrix or a positive scalar, d is a 1x2 vector )
% with three boundary conditions: 
% 1. Dirichlet: u = gD
% 2. Neumann: c*grad(u)·n = gN
% 3. Robin: c*grad(u)·n + gR1*u = gR2
% Remark: the mesh data is got from iFEM software package


% setpath;
close all;
clear;
% clc;
tic;

exnum = 4; % example number
elemtype = 3; % 3 for triangle mesh and 4 for quadrilateral mesh
femtype.type = 'pgm'; %PG method
%femtype.type = 'ipdg'; %IPDG method

trbasis_type = 'l'; trbasis_order = 2;
tsbasis_type = 'l'; tsbasis_order = 2;
Gauss_order = 7;
errtype = 'l2';
plottype = '1'; % '0':draw mesh; '1':plot error; '2':plot solution;
                  % '3':stiffness matrix; '4':error surface; '5':cpu time
plotdx_order = 0;
plotdy_order = 0;
upperbd = 1e2;  
refine = 8;

switch errtype
    case 'linf'
        ord = trbasis_order + 1; soblev = 'L^{\infty}';
    case 'l2'
        ord = trbasis_order + 1; soblev = 'L^2'; 
    case 'h1'
        ord = trbasis_order;  soblev = 'H^1'; 
    case 'h1semi'
        ord = trbasis_order;  soblev = 'H^1 semi';
    case 'h2'
        ord = trbasis_order - 1;  soblev = 'H^2';
    case 'h2semi'
        ord = trbasis_order - 1;  soblev = 'H^2 semi';
    case 'h3'
        ord = trbasis_order - 2;  soblev = 'H^3';
    case 'h3semi'
        ord = trbasis_order - 2;  soblev = 'H^3 semi';
    case 'h4'
        ord = trbasis_order - 3;  soblev = 'H^4';
    case 'h4semi'
        ord = trbasis_order - 3;  soblev = 'H^4 semi';
    case 'h5'
        ord = trbasis_order - 4;  soblev = 'H^5';
    case 'h5semi'
        ord = trbasis_order - 4;  soblev = 'H^5 semi';
    case 'h6'
        ord = trbasis_order - 5;  soblev = 'H^6';
    case 'h6semi'
        ord = trbasis_order - 5;  soblev = 'H^6 semi';
end


switch exnum
    %% example -1: a polynomial case, \Omega=(-1,1)x(-1,1)
    case -1    
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2); % structure mesh
                %[node,elem] = regpolygon(5,1);  % polygon mesh
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = 1;
        coeff.d = [1,1];
        coeff.k = -1;
        f = @(x,y) -4 + 3*x + 3*y - (x.^2 + x.*y + y.^2);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) x.^2 + x.*y + y.^2; % exact solution
        uclass{2,1} = @(x,y) x + 2*y; % uy
        uclass{2,2} = @(x,y) 2*x + y; % ux
        uclass{3,1} = @(x,y) 2; % uyy
        uclass{3,2} = @(x,y) 1; % uxy
        uclass{3,3} = @(x,y) 2; % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 0: a polynomial case, \Omega=(-1,1)x(-1,1)
    case 0    
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2); % structure mesh
                %[node,elem] = regpolygon(5,1);  % polygon mesh
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = 0;
        f = @(x,y) -20*x.^3.*y-2;
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) x.^5.*y + y.^2; % exact solution
        uclass{2,1} = @(x,y) x.^5 + 2*y; % uy
        uclass{2,2} = @(x,y) 5*x.^4.*y; % ux
        uclass{3,1} = @(x,y) 2+x-x; % uyy
        uclass{3,2} = @(x,y) 5*x.^4; % uxy
        uclass{3,3} = @(x,y) 20*x.^3.*y; % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 1: a simple case for test, \Omega=(-1,1)x(-1,1)
    case 1    
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2); % structure mesh
%                 g = [ 2     2     2     2
%                     -1     1     1    -1
%                     1     1    -1    -1
%                     1     1    -1    -1
%                     1    -1    -1     1
%                     0     0     0     0
%                     1     1     1     1 ];
%                 [p,e,t] = initmesh(g,'hmax',2);
%                 node = p'; elem = t(1:3,:)';  % non-structure mesh
                %[node,elem] = regpolygon(5,1);  % polygon mesh
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = 0;
        f = @(x,y) -2*exp(x+y);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) exp(x+y); % exact solution
        uclass{2,1} = @(x,y) exp(x+y); % uy
        uclass{2,2} = @(x,y) exp(x+y); % ux
        uclass{3,1} = @(x,y) exp(x+y); % uyy
        uclass{3,2} = @(x,y) exp(x+y); % uxy
        uclass{3,3} = @(x,y) exp(x+y); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};
        % 边界必须不相交,也可以写为
%         bdinfo.function{1} = @(x,y) exp(-1+y).*(x==-1 & y~=1) + ...
%                             exp(1+y).*(x==1 & y~=-1) + ...
%                             exp(x-1).*(y==-1 & x~=-1) + ...
%                             exp(x+1).*(y==1 & x~=1);

    %% example 2: Poisson equation, \Omega=(-1,1)x(-1,1)
    case 2    
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2); % structure mesh
                %load lakemesh.mat; node = node/5; % lake mesh 
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = 0;
        f = @(x,y) -y.*(1-y).*(1-x-x.*x/2).*exp(x+y)-x.*(1-x/2).*(-3*y-y.*y).*exp(x+y);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) x.*y.*(1-x/2).*(1-y).*exp(x+y); % exact solution
        uclass{2,1} = @(x,y) 1/2*exp(x + y).*(-2 + x).*x.*(-1 + y + y.^2); % uy
        uclass{2,2} = @(x,y) 1/2*exp(x + y).*(-2 + x.^2).*(-1 + y).*y;   % ux
        uclass{3,1} = @(x,y) 1/2*exp(x + y).*(-2 + x).*x.*y.*(3 + y); % uyy
        uclass{3,2} = @(x,y) 1/2*exp(x + y).*(-2 + x.^2).*(-1 + y + y.^2); % uxy
        uclass{3,3} = @(x,y) 1/2*exp(x + y).*(-2 + 2*x + x.^2).*(-1 + y).*y; % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 3: non-symmetric elliptic equation, \Omega=(-1,1)x(-1,1)
    case 3   
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2);
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = 1;
        coeff.d = cell(1,2); coeff.d{1} = @(x,y) x; coeff.d{2} = @(x,y) y;
        coeff.k = @(x,y) -(x+y);
        f = @(x,y) -2*exp(x+y);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) exp(x+y); % exact solution
        uclass{2,1} = @(x,y) exp(x+y); % uy
        uclass{2,2} = @(x,y) exp(x+y); % ux
        uclass{3,1} = @(x,y) exp(x+y); % uyy
        uclass{3,2} = @(x,y) exp(x+y); % uxy
        uclass{3,3} = @(x,y) exp(x+y); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 4: non-isotropic elliptic equation
    case 4    
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2);
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = cell(2,2);
        coeff.c{1,1} = @(x,y) exp(x);  coeff.c{1,2} = @(x,y) -exp(2*x);
        coeff.c{2,1} = @(x,y) exp(2*y);  coeff.c{2,2} = @(x,y) exp(y);
        coeff.d = cell(1,2); 
        coeff.d{1} = @(x,y) exp(x)+2*exp(2*y); 
        coeff.d{2} = @(x,y) -2*exp(2*x)+exp(y);
        coeff.k = @(x,y) -exp(x)-exp(y);
        f = @(x,y) (-exp(2*x)+exp(2*y)).*sin(y).*cos(x);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) sin(x).*cos(y); % exact solution
        uclass{2,1} = @(x,y) -sin(x).*sin(y); % uy
        uclass{2,2} = @(x,y) cos(x).*cos(y); % ux
        uclass{3,1} = @(x,y) -sin(x).*cos(y); % uyy
        uclass{3,2} = @(x,y) -cos(x).*sin(y); % uxy
        uclass{3,3} = @(x,y) -sin(x).*cos(y); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 5: multiscale elliptic equation, \Omega=(-1,1)x(-1,1)
    case 5    
        e = 0.01;
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2);
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = @(x,y) cos(x.*y/e)+1;
        coeff.d = [0,0];
        coeff.k = 0;
        f = @(x,y) (x.^2+y.^2)/e^2.*(sin(2*x.*y/e)+sin(x.*y/e));
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) sin(x.*y/e); % exact solution
        uclass{2,1} = @(x,y) x/e.*cos(x.*y/e); % uy
        uclass{2,2} = @(x,y) y/e.*cos(x.*y/e); % ux
        uclass{3,1} = @(x,y) -x.^2/e^2.*sin(x.*y/e); % uyy
        uclass{3,2} = @(x,y) cos(x.*y/e)/e - x.*y.*sin(x.*y/e)/e^2; % uxy
        uclass{3,3} = @(x,y) -y.^2/e^2.*sin(x.*y/e); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 6: another multiscale elliptic equation, \Omega=(-1,1)x(-1,1)
    case 6    
        e1 = 10*pi; e2 = 100*pi;
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],2);
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],2);
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = 0;
        f = @(x,y) 4*e1^2*x.^2.*sin(e1*x.^2) - 2*e1*cos(e1*x.^2) + ...
            4*e2^2*y.^2.*sin(e2*y.^2) - 2*e2*cos(e2*y.^2);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) sin(e1*x.^2)+sin(e2*y.^2); % exact solution
        uclass{2,1} = @(x,y) 2*e2*y.*cos(e2*y.^2); % uy
        uclass{2,2} = @(x,y) 2*e1*x.*cos(e1*x.^2); % ux
        uclass{3,1} = @(x,y) 2*e2*cos(e2*y.^2)-4*e2^2*y.^2.*sin(e2*y.^2); % uyy
        uclass{3,2} = @(x,y) 0+x-x; % uxy
        uclass{3,3} = @(x,y) 2*e1*cos(e1*x.^2)-4*e1^2*x.^2.*sin(e1*x.^2); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};

    %% example 7: L-shape, \Omega=(-1,1)x(-1,1) \ (0,1)x(-1,0)
    case 7   
        switch elemtype
            case 3
                [node,elem] = squaremesh([-1,1,-1,1],1);
                [node,elem] = delmesh(node,elem,'x>0 & y<0');
            case 4
                [node,elem] = squarequadmesh([-1,1,-1,1],1);
                [node,elem] = delmesh(node,elem,'x>0 & y<0');
        end
        mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
        bdFlag = setboundary(node,elem,'Dirichlet');
        bdinfo = genMeshbdinfo(elem,bdFlag);
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = @lshape;
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = 0;
        f = 0;
        uclass = cell(7,9);
        uclass{1,1} = @lshape; % exact solution
        uclass{2,1} = @lshapedy; % uy
        uclass{2,2} = @lshapedx; % ux
        uclass{3,1} = @lshapeddyy; % uyy
        uclass{3,2} = @lshapeddxy; % uxy
        uclass{3,3} = @lshapeddxx; % uxx  

    %% example 8: complex Helmholtz equation I: plane wave, \Omega=(-1,1)x(-1,1)
    case 8    
        k = 50;  theta = pi;
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
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = -k^2;
        f = 0;
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) exp(1i*k*(x*cos(theta)+y*sin(theta))); % exact solution
        uclass{2,1} = @(x,y) 1i*k*sin(theta).*exp(1i*k*(x*cos(theta)+y*sin(theta))); % uy
        uclass{2,2} = @(x,y) 1i*k*cos(theta).*exp(1i*k*(x*cos(theta)+y*sin(theta))); % ux
        uclass{3,1} = @(x,y) -k^2*sin(theta)^2.*exp(1i*k*(x*cos(theta)+y*sin(theta))); % uyy
        uclass{3,2} = @(x,y) -k^2*sin(theta)*cos(theta).*exp(1i*k*(x*cos(theta)+y*sin(theta))); % uxy
        uclass{3,3} = @(x,y) -k^2*cos(theta)^2.*exp(1i*k*(x*cos(theta)+y*sin(theta))); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1};  

    %% example 9: complex Helmholtz equation II, \Omega=(-1,1)x(-1,1)
    case 9    
        k = 100; 
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
        coeff.c = 1;
        coeff.d = [0,0];
        coeff.k = -k^2;
        f = @(x,y) sin(k*sqrt(x.^2+y.^2))./sqrt(x.^2+y.^2);
        uclass = cell(7,9);
        uclass{1,1} = @(x,y) cos(k*sqrt(x.^2+y.^2))/k-(cos(k)+1i*sin(k))*besselj(0,k*sqrt(x.^2+y.^2))...
                             /k/(besselj(0,k)+1i*besselj(1,k)); % exact solution
        uclass{2,1} = @(x,y) ( -sin(k*sqrt(x.^2+y.^2))+(cos(k)+1i*sin(k))*besselj(1,k*sqrt(x.^2+y.^2))...
                             /(besselj(0,k)+1i*besselj(1,k)) ).*y./sqrt(x.^2+y.^2); % uy
        uclass{2,2} = @(x,y) ( -sin(k*sqrt(x.^2+y.^2))+(cos(k)+1i*sin(k))*besselj(1,k*sqrt(x.^2+y.^2))...
                             /(besselj(0,k)+1i*besselj(1,k)) ).*x./sqrt(x.^2+y.^2); % ux
        uclass{3,1} = @(x,y) ( -sin(k*sqrt(x.^2+y.^2))+(cos(k)+1i*sin(k))*besselj(1,k*sqrt(x.^2+y.^2))...
                             /(besselj(0,k)+1i*besselj(1,k)) ).*(1./sqrt(x.^2+y.^2)-y.^2./sqrt(x.^2+y.^2).^3) ...
                             +(y.^2./sqrt(x.^2+y.^2).^2).*(-k*cos(k*sqrt(x.^2+y.^2))+k*(cos(k)+1i*sin(k))/2/(besselj(0,k)+1i*besselj(1,k)) ...
                             *(besselj(0,k*sqrt(x.^2+y.^2))-besselj(2,k*sqrt(x.^2+y.^2))) ); % uyy
        uclass{3,2} = @(x,y) (x.*y./sqrt(x.^2+y.^2).^2).*(-k*cos(k*sqrt(x.^2+y.^2))+k*(cos(k)+1i*sin(k))/2/(besselj(0,k)+1i*besselj(1,k)) ...
                             *(besselj(0,k*sqrt(x.^2+y.^2))-besselj(2,k*sqrt(x.^2+y.^2))) ); % uxy
        uclass{3,3} = @(x,y) ( -sin(k*sqrt(x.^2+y.^2))+(cos(k)+1i*sin(k))*besselj(1,k*sqrt(x.^2+y.^2))...
                             /(besselj(0,k)+1i*besselj(1,k)) ).*(1./sqrt(x.^2+y.^2)-x.^2./sqrt(x.^2+y.^2).^3) ...
                             +(x.^2./sqrt(x.^2+y.^2).^2).*(-k*cos(k*sqrt(x.^2+y.^2))+k*(cos(k)+1i*sin(k))/2/(besselj(0,k)+1i*besselj(1,k)) ...
                             *(besselj(0,k*sqrt(x.^2+y.^2))-besselj(2,k*sqrt(x.^2+y.^2))) ); % uxx
        bdinfo.function = cell(1,4);
        bdinfo.function{1} = uclass{1,1}; 

end

Efem = zeros(1,refine);  % FEM error: ||u-u_h||
EI = zeros(1,refine); % interpolation error: ||u-I_h(u)||
EIfem = zeros(1,refine); % ||I_h(u)-u_h||
T = zeros(1,refine);  % cpu time
H = zeros(1,refine);  % mesh size
Dof = zeros(1,refine);  % degree of freedom
uzeros = cell(7,9);
for i = 1:7
    for j = 1:9
        uzeros{i,j} = @(x,y) 0+x-x;
    end
end

for i = 1:refine
    t1=cputime;
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
    
    % compute mesh size
    e1 = mesh.P(:,mesh.T(1,:))-mesh.P(:,mesh.T(2,:));
    e2 = mesh.P(:,mesh.T(2,:))-mesh.P(:,mesh.T(3,:));
    e3 = mesh.P(:,mesh.T(3,:))-mesh.P(:,mesh.T(1,:));
    h0 = min([e1(1,:).^2+e1(2,:).^2, e2(1,:).^2+e2(2,:).^2, e3(1,:).^2+e3(2,:).^2]);
    H(i) = sqrt(h0);

    %% FEM main routine
    [uh,~,basis,A,F] = lineEllip2(mesh,coeff,f,bdinfo,trbasis_type,...
        trbasis_order,tsbasis_type,tsbasis_order,Gauss_order,femtype);
    % FEM interpolation
    uI = uclass{1,1}(basis.trial.Pb(1,:)',basis.trial.Pb(2,:)');

    switch elemtype
        case 3
            Efem(i) = triError2(mesh,basis,uclass,uh,errtype,Gauss_order,femtype);
            EI(i) = triError2(mesh,basis,uclass,uI,errtype,Gauss_order,femtype);
            EIfem(i) = triError2(mesh,basis,uzeros,uI-uh,errtype,Gauss_order,femtype);
        case 4
            Efem(i) = quadError2(mesh,basis,uclass,uh,errtype,Gauss_order,femtype);
            EI(i) = quadError2(mesh,basis,uclass,uI,errtype,Gauss_order,femtype);
            EIfem(i) = triError2(mesh,basis,uzeros,uI-uh,errtype,Gauss_order,femtype);
    end

    t2=cputime;
    T(i) = t2-t1;
    Dof(i) = basis.trial.Nb;
    fprintf('i = %d, dof = %d \n',i,Dof(i));
end

toc;

if contains(plottype,'0') % draw mesh
    figure(10);
    showmesh(node,elem);
    %findnode(basis.test.Pb');
    %findelem(node,elem);
end

if contains(plottype,'1') % plot error
    R = Efem;
    switch isnan(R(1))
        case 1
            R(1) = 0;
        case 0
            R(1) = 0.5*R(1);
    end
    for i=2:length(Efem)
        R(i) = R(i-1)/2^ord;
    end
    figure(1); 
%     subplot(121);
%     set(gcf,'unit','normalized','position',[0.2, 0.2, 0.8, 0.6]);
%     semilogy(log(1./H),Efem,'-*',log(1./H),EI,'-o',log(1./H),EIfem,'-.x',...
%              log(1./H),R,'--+','linewidth',1); hold on;
%     str1 = ['Lagrange ',num2str(trbasis_order),' order element ',soblev,' error'];
%     str2 = strcat(['$||u-u_h||_{',soblev,'}$']);
%     str3 = strcat(['$||u-I_h(u)||_{',soblev,'}$']);
%     str4 = strcat(['$||I_h(u)-u_h||_{',soblev,'}$']);
%     str5 = ['$y=O(h^',num2str(ord),')$'];
%     title(str1);
%     legend(str2,str3,str4,str5,'location','southwest','interpreter','latex');
%     xlabel('$-\log h$','Interpreter','latex');
%     ylabel('Error','Interpreter','latex');
%     set(gca,'fontsize',18);
% 
%     subplot(122);
    set(gcf,'unit','normalized','position',[0.2, 0.2, 0.5, 0.7]);
    loglog(Dof,Efem,'-*',Dof,EI,'-o',Dof,EIfem,'-.x',Dof,R,'--+','linewidth',1); hold on;
    str1 = ['Lagrange ',num2str(trbasis_order),' order element ',soblev,' error'];
    str2 = strcat(['$||u-u_h||_{',soblev,'}$']);
    str3 = strcat(['$||u-I_h(u)||_{',soblev,'}$']);
    str4 = strcat(['$||I_h(u)-u_h||_{',soblev,'}$']);
    str5 = ['$y=O(h^',num2str(ord),')$'];
    title(str1);
    legend(str2,str3,str4,str5,'location','southwest','interpreter','latex');
    xlabel('Dof','Interpreter','latex');
    ylabel('Error','Interpreter','latex');
    set(gca,'fontsize',18);
end

if contains(plottype,'2') % plot fem and exact solution
    s = num2str(plotdx_order); t = num2str(plotdy_order);
    st = num2str(plotdx_order + plotdy_order);
    str1 = strcat(['FEM linear interpolation solution $I_h(d^',st,'u_h',...
                    '/dx^',s,'dy^',t,')$']);
    str2 = strcat(['Exact solution $I_h(d^',st,'u/dx^',s,'dy^',t,')$']);
    figure(2);
    subplot(121); % plot fem solution 
    set(gcf,'unit','normalized','position',[0.2, 0.2, 0.9, 0.6]);
%     (only works when \Omega is convex polygon)
%     a = 30; b = 30; 
%     xmin = min(mesh.P(1,:));  xmax = max(mesh.P(1,:));
%     ymin = min(mesh.P(2,:));  ymax = max(mesh.P(2,:));
%     [X,Y] = meshgrid(xmin:(xmax-xmin)/a:xmax, ymin:(ymax-ymin)/a:ymax);
%     X = reshape(X,1,[]); Y = reshape(Y,1,[]); bdnode = bdinfo.nodes(2,:);
%     plotpt = [X,mesh.P(1,bdnode); Y,mesh.P(2,bdnode)]; % include boundary modes
%     [Iuh, outindex] = utriFEM2(mesh,basis,uh,plotpt(1,:)',plotpt(2,:)',s,t,femtype,upperbd);
%     Iuh(outindex) = []; plotpt(:,outindex) = [];
%     tri = delaunay(plotpt(1,:),plotpt(2,:));
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
        %trisurf(tri,plotpt(1,:),plotpt(2,:),Iuh); shading interp;
        pdesurf(p0,t0,Iuh);
    else
        %trisurf(tri,plotpt(1,:),plotpt(2,:),abs(Iuh)); shading interp;
        %trisurf(tri,plotpt(1,:),plotpt(2,:),-angle(Iuh)); shading interp;
        pdesurf(p0,t0,abs(Iuh));
        %pdesurf(p0,t0,-angle(Iuh));
    end
    view(0,90);
    colorbar;
    title(str1,'Interpreter','latex');
    set(gca,'fontsize',18);
    grid on;

    subplot(122); % plot exact solution (linear interpolation)
    uexact = uclass{plotdx_order + plotdy_order + 1, plotdx_order + 1}...
            (mesh.P(1,:)',mesh.P(2,:)');
    p0 = mesh.P;
    t0 = [mesh.T;ones(1,length(mesh.T(1,:)))];
    if isreal(uexact)
        pdesurf(p0,t0,uexact);
    else
        pdesurf(p0,t0,abs(uexact));
        %pdesurf(p0,t0,angle(uexact));
    end
    shading interp;
    colormap parula;
    view(0,90);
    colorbar;
    title(str2,'Interpreter','latex');
    set(gca,'fontsize',18);
    grid on;
end

if contains(plottype,'3') % plot stiffness matrix
    figure(3)
    spy(A);
    title("Nonzero element in stiffness matrix");
    set(gca,'fontsize',18);
end

if contains(plottype,'4') % plot absolute error surface
    figure(4)
    pdesurf(p0,t0,abs(uI - uh));
    colormap parula;
    view(0,90);
    colorbar;
    title('$||I_h(u)-u_h||$','Interpreter','latex');
    set(gca,'fontsize',18);
    grid on;
end

if contains(plottype,'5') % plot cpu time during each refinement
    figure(5)
    plot(1:1:length(T),T,'-o','linewidth',1);
    title("CPU time");
    set(gca,'fontsize',18);
end

 