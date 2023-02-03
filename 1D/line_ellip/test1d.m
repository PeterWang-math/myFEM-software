% This script is a test for FEM solver of 1d stable elliptic equation:
%       -(c*u')'+d*u'+k*u=f on \Omega, 
% ( c>0 is a positive scalar)
% with three boundary conditions: 
% 1. Dirichlet: u=g
% 2. Neumann: c*u'=r
% 3. Robin: c*u'+p·u=q

tic;

exnum = 5;  % the example number, Don't forget to change this ! ! ! ! ! !
femtype.type = 'pgm'; %PG method
%femtype.type = 'ipdg'; %IPDG method
mesh_type=1;
trbasis_type='l'; trbasis_order=1;
tsbasis_type='l'; tsbasis_order=1;
Gauss_nodes = 7;
errtype = 'linf';
refine = 10;

bdinfo.nodes = zeros(2,2);

switch exnum
    case 1    % example 1: symmetric elliptic equation, \Omega=(0,1)
        a=0; b=1;
        coeff.c = @(x) exp(x);
        coeff.d = @(x) 0;
        coeff.k = @(x) 0;
        f = @(x) -exp(x).*(cos(x)-2*sin(x)-x.*cos(x)-x.*sin(x));
        uclass = cell(7,1);
        bdinfo.nodes(1,:) = [1, 1];  % set Dirichlet boundary condition type
        bdinfo.value = [0, cos(1)];
        uclass{1} = @(x) cos(x).*x;  % exact solution

        % pure Neumann b.c.
        % bdinfo.nodes(1,:) = [2, 2];
        % bdinfo.value = [1, cos(1)-sin(1)];
        % uclass{1} = @(x) cos(x).*x+1-sin(1)-cos(1);

        uclass{2} = @(x) cos(x)-sin(x).*x;
        uclass{3} = @(x) -2*sin(x)-cos(x).*x;
        uclass{4} = @(x) -3*cos(x)+sin(x).*x;
        uclass{5} = @(x) 4*sin(x)+cos(x).*x;
        uclass{6} = @(x) 5*cos(x)-sin(x).*x;
        uclass{7} = @(x) -6*sin(x)-cos(x).*x;

    case 2    % example 2: non-symmetric elliptic equation, \Omega=(0,1)
        a=0; b=1;
        coeff.c = @(x) exp(-x);
        coeff.d = @(x) 1;
        coeff.k = @(x) 1;
        f = @(x) (2*x+1).*exp(x)-1;
        bdinfo.nodes(1,:) = [3, 1]; % set Robin & Dirichlet boundary condition type
        bdinfo.value = [1, exp(1); 1, 3*exp(1)];
        uclass = cell(7,1);
        uclass{1} = @(x) exp(x).*x;  % exact solution
        uclass{2} = @(x) exp(x).*(x+1);
        uclass{3} = @(x) exp(x).*(x+2);
        uclass{4} = @(x) exp(x).*(x+3);
        uclass{5} = @(x) exp(x).*(x+4);
        uclass{6} = @(x) exp(x).*(x+5);
        uclass{7} = @(x) exp(x).*(x+6);

    case 3    % example 3: multiscale elliptic equation, \Omega=(0,pi)
        e = 0.05; a=0; b=pi;
        coeff.c = @(x) 1;
        coeff.d = @(x) (2/e*cos(x).*sin(x/e)+sin(x).*cos(x/e))./(cos(x).*cos(x/e));
        coeff.k = @(x) 0;
        f = @(x) -exp(-x).*cos(x/e).*(-2*e + (-3 + e^2)*tan(x/e) + 2*e*tan(x/e).^2 + ...
            e*tan(x).*(-1 + e*tan(x/e)))/e^2;
        bdinfo.nodes(1,:) = [1, 1];  % set Dirichlet boundary condition type
        bdinfo.value = [0, exp(-pi)*sin(pi/e)];
        uclass = cell(7,1);
        uclass{1} = @(x) exp(-x).*sin(x/e);  % exact solution
        uclass{2} = @(x) exp(-x).*(cos(x/e)/e - sin(x/e));
        uclass{3} = @(x) exp(-x).*(-2*e*cos(x/e) + (-1 + e^2)*sin(x/e))/e^2;
        uclass{4} = @(x) exp(-x).*((-1 + 3*e^2)*cos(x/e) - e*(-3 + e^2)*sin(x/e))/e^3;
        uclass{5} = @(x) exp(-x).*(-4*e*(-1 + e^2)*cos(x/e) + (1 - 6*e^2 + e^4)*sin(x/e))/e^4;
        uclass{6} = @(x) exp(-x).*((1 - 10*e^2 + 5*e^4)*cos(x/e) - ...
            e*(5 - 10*e^2 + e^4)*sin(x/e))/e^5;
        uclass{7} = @(x) exp(-x).*(-2*e*(3 - 10*e^2 + 3*e^4)*cos(x/e) + ...
            (-1 + 15*e^2 - 15*e^4 + e^6)*sin(x/e))/e^6;

    case 4    % example 4: another multiscale elliptic equation, \Omega=(-3,3)
        e = 500; a=0; b=3;
        coeff.c = @(x) 1;
        coeff.d = @(x) 0;
        coeff.k = @(x) 0;
        bdinfo.nodes(1,:) = [1, 1];  % set Dirichlet boundary condition type
        bdinfo.value = [0, 0]; 
        uclass = cell(7,1);
        uclass{1} = @(x) sin(2*pi*x).*(2+0.1*cos(e*pi*x)).*exp(-x.^2);  % exact solution
        uclass{2} = @(x) 2*pi*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         2*x.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         (e*pi*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/10;
        uclass{3} = @(x) 4*x.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         2*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         4*pi^2*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         (2*e*pi^2*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 - ...
                         (e^2*pi^2*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/10 - ...
                         8*x*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         (2*e*x*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5;
        uclass{4} = @(x) 12*x.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 12*pi*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         8*x.^3.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 8*pi.^3*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         (e.^3*pi.^3*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/10 + (3*e*pi*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 + ...
                         24*x.^2*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + 24*x*pi.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         (6*e*pi.^3*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 - (3*e.^2*pi.^3*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x))/5 + ...
                         (12*e*x*pi.^2.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 - (6*e*x.^2*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 + ...
                         (3*e.^2*x*pi.^2.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/5;
        uclass{5} = @(x) 12*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 48*x.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         16*x.^4.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + 48*pi.^2*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         16*pi.^4*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + 64*x*pi.^3.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         64*x.^3*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + (24*e*pi.^2*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 + ...
                         (16*e*pi.^4*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 - 96*x.^2*pi.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         (6*e.^2*pi.^2*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/5 + (12*e.^2*pi.^4*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/5 + ...
                         (4*e.^3*pi.^4*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 + (e.^4*pi.^4*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/10 + ...
                         96*x*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - (48*e*x*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 + ...
                         (16*e*x.^3*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 + (24*e.^2*x*pi.^3.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x))/5 - ...
                         (48*e*x.^2*pi.^2.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 - (4*e.^3*x*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 - ...
                         (24*e*x*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 - (12*e.^2*x.^2*pi.^2.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/5;
        uclass{6} = @(x) 160*x.^3.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 32*x.^5.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         120*pi*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - 120*x.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         160*pi.^3*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + 32*pi.^5*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         2*e.^3*pi.^3*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - 4*e.^3*pi.^5*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - ...
                         (e.^5*pi.^5*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/10 - 6*e*pi*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - ...
                         480*x.^2*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + 160*x.^4*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         480*x*pi.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 160*x*pi.^4.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         24*e*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - 8*e*pi.^5*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - ...
                         320*x.^2*pi.^3.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + 320*x.^3*pi.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         12*e.^2*pi.^3.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) + 8*e.^2*pi.^5.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) + ...
                         e.^4*pi.^5*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) - 48*e*x*pi.^2.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - ...
                         32*e*x*pi.^4.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) + 24*e*x.^2*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - ...
                         8*e*x.^4*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - 12*e.^2*x*pi.^2.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) + ...
                         32*e*x.^3*pi.^2.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - 24*e.^2*x*pi.^4.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - ...
                         8*e.^3*x*pi.^4.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - e.^4*x*pi.^4.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) + ...
                         48*e*x.^2*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - 24*e.^2*x.^2*pi.^3.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) + ...
                         8*e.^2*x.^3*pi.^2.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) + 4*e.^3*x.^2*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x);
        uclass{7} = @(x) 720*x.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 120*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         480*x.^4.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + 64*x.^6.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         720*pi.^2*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 480*pi.^4*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         64*pi.^6*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 1920*x*pi.^3.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         1920*x.^3*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - 384*x*pi.^5.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         384*x.^5*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) - 72*e*pi.^2*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - ...
                         96*e*pi.^4*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - (96*e*pi.^6*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 + ...
                         1280*x.^3*pi.^3.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + 2880*x.^2*pi.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         960*x.^2*pi.^4.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - 960*x.^4*pi.^2.*exp(-x.^2).*sin(2*pi*x).*(cos(pi*e*x)/10 + 2) - ...
                         18*e.^2*pi.^2*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - 72*e.^2*pi.^4*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - ...
                         24*e.^3*pi.^4*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - 24*e.^2*pi.^6*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - ...
                         3*e.^4*pi.^4*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - 16*e.^3*pi.^6*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - ...
                         6*e.^4*pi.^6*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - (6*e.^5*pi.^6*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x))/5 - ...
                         (e.^6*pi.^6*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x))/10 - 1440*x*pi.*exp(-x.^2).*cos(2*pi*x).*(cos(pi*e*x)/10 + 2) + ...
                         288*e*x*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) - 96*e*x.^3*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) + ...
                         96*e*x*pi.^5.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) + (96*e*x.^5*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 - ...
                         144*e.^2*x*pi.^3.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) - 96*e.^2*x*pi.^5.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) - ...
                         12*e.^4*x*pi.^5.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) + 288*e*x.^2*pi.^2.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) + ...
                         192*e*x.^2*pi.^4.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - 96*e*x.^4*pi.^2.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) - ...
                         192*e*x.^3*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) + 24*e.^3*x*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) + ...
                         48*e.^3*x*pi.^5.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) + (6*e.^5*x*pi.^5.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x))/5 + ...
                         96*e.^2*x.^3*pi.^3.*exp(-x.^2).*cos(2*pi*x).*cos(pi*e*x) + 72*e*x*pi.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x) + ...
                         72*e.^2*x.^2*pi.^2.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) + 144*e.^2*x.^2*pi.^4.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - ...
                         24*e.^2*x.^4*pi.^2.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) + 48*e.^3*x.^2*pi.^4.*sin(pi*e*x).*exp(-x.^2).*cos(2*pi*x) + ...
                         6*e.^4*x.^2*pi.^4.*exp(-x.^2).*sin(2*pi*x).*cos(pi*e*x) - 16*e.^3*x.^3*pi.^3.*sin(pi*e*x).*exp(-x.^2).*sin(2*pi*x);
        f = @(x) -uclass{3}(x);

    case 5    % example 5: discontinues coefficient, \Omega=(-1,1)
        a=-1; b=1; eps=1e-2;
        %coeff.c = @(x) 1*(x<-eps)+2*(x>eps)+...
        %               (1/(2*eps)*x+3/2).*((x<=eps)&(x>-eps));
        coeff.c = @(x) 1*(x<0)+2*(x>=0);
        coeff.d = @(x) 0;
        coeff.k = @(x) 0;
        %f = @(x) -exp(x).*(x<-eps)-2*exp(x).*(x>eps)...
        %         -exp(x).*(1/(2*eps)+3/2+1/(2*eps)*x).*((x<=eps)&(x>-eps));
        f = @(x) -exp(x).*(x<0)-2*exp(x).*(x>=0);
        bdinfo.nodes(1,:) = [1, 1];  % set Dirichlet boundary condition type
        bdinfo.value = [exp(-1), exp(1)];
        uclass = cell(7,1);
        uclass{1} = @(x) exp(x);  % exact solution
        uclass{2} = @(x) exp(x);
        uclass{3} = @(x) exp(x);
        uclass{4} = @(x) exp(x);
        uclass{5} = @(x) exp(x);
        uclass{6} = @(x) exp(x);
        uclass{7} = @(x) exp(x);
        
end

E = zeros(1,refine); E1=E; E2=E; E3=E; E4=E; E5=E; E6=E; E7=E;
for i = 1:refine
    mesh_N=2^i;
    mesh = genmesh1(a,b,mesh_type,mesh_N);
    femtype.ipdg.gamma0 = 1*ones(1,length(mesh.P));
    femtype.ipdg.gamma1 = 1*ones(1,length(mesh.P));
    femtype.ipdg.alpha = 2*ones(1,length(mesh.P));
    femtype.ipdg.beta = -1*ones(1,length(mesh.P));
    femtype.ipdg.theta = 1; % SIPG
    %femtype.ipdg.theta = -1; % NIPG
    %femtype.ipdg.theta = 0; % IIPG
    bdinfo.nodes(2,:) = [1,mesh_N+1];
    [uh,ufem,basis,A,F] = line_ellip1(mesh,coeff,f,bdinfo,...
        trbasis_type,trbasis_order,tsbasis_type,tsbasis_order,...
        Gauss_nodes,femtype);
    E(i) = get_error1(mesh,basis,uclass,uh,errtype,Gauss_nodes,femtype);
    E1(i) = get_error1(mesh,basis,uclass,uh,'l2',Gauss_nodes,femtype);
    E2(i) = get_error1(mesh,basis,uclass,uh,'h1',Gauss_nodes,femtype);
    E3(i) = get_error1(mesh,basis,uclass,uh,'h2',Gauss_nodes,femtype);
    E4(i) = get_error1(mesh,basis,uclass,uh,'h3',Gauss_nodes,femtype);
    E5(i) = get_error1(mesh,basis,uclass,uh,'h4',Gauss_nodes,femtype);
    E6(i) = get_error1(mesh,basis,uclass,uh,'h5',Gauss_nodes,femtype);
    E7(i) = get_error1(mesh,basis,uclass,uh,'h6',Gauss_nodes,femtype);
    fprintf('i=%d, dof=%d \n',i,basis.trial.Nb);
end

toc;

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

switch trbasis_type
    case 'l'
        basis_type = 'Lagrange ';
    case 'h'
        basis_type = 'Hermite ';
end

figure(1); % plot error
R = E;
switch isnan(R(1))
    case 1
        R(1) = 0;
    case 0
        R(1) = R(1) + 0.01;
end
for i=2:length(E)
    R(i) = R(i-1)/2^ord;
end
semilogy(1:1:length(E),E,'-o',1:1:length(R),R,'--','linewidth',1); hold on;
str1 = [basis_type,num2str(trbasis_order),' order element ',soblev,' error'];
str2 = ['y=O(h^',num2str(ord),')'];
legend(str1,str2,'location','southwest','interpreter','tex');
set(gca,'fontsize',18);

figure(2); % plot exact solution and fem solution
x=linspace(a,b,1e4); uf=x; u=x;
for i=1:length(x)
    uf(i) = ufem(x(i));
    u(i) = uclass{1}(x(i));
end
plot(x,uf,'-o','linewidth',1); grid on; hold on;
plot(x,u,'linewidth',1);
legend([basis_type,num2str(trbasis_order),' order element'],...
    'Exact solution','location','southwest');
set(gca,'fontsize',18);
 

 