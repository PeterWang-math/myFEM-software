% a test script for solving eigenvalue problem of Laplace equaction: 
%           -\Delta u = \lambda u
%           u = 0 on \gamma
%           ||u||_{L^2} = 1

close all;
%clear;
% clc;
tic;

k = 100; % return the smallest k eigenvalue
elemtype = 3; % 3 for triangle mesh and 4 for quadrilateral mesh
femtype.type = 'pgm'; %PG method
trbasis_type = 'l'; trbasis_order = 1;
tsbasis_type = 'l'; tsbasis_order = 1;
Gauss_order = 7;
errtype = 'linf';
plottype = ''; % '0':draw mesh; '1':plot error; '2':plot solution;
                  % '3':stiffness matrix; '4':error surface; '5':cpu time
plotdx_order = 0;
plotdy_order = 0;
upperbd = 1e2;  


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

% Initial mesh and bdinfo
switch elemtype
    case 3
        %[node,elem] = squaremesh([-1,1,-1,1],2); % structure mesh
%         g = [ 2     2     2     2
%             -1     1     1    -1
%             1     1    -1    -1
%             1     1    -1    -1
%             1    -1    -1     1
%             0     0     0     0
%             1     1     1     1 ];
%         [p,e,t] = initmesh(g); 
        node = p'; elem = t(1:3,:)';  % non-structure mesh
        %[node,elem] = regpolygon(5,1);  % polygon mesh
    case 4
        [node,elem] = squarequadmesh([-1,1,-1,1],2);
end
mesh.P = node';  mesh.T = elem';  mesh.N = size(elem,1);
bdFlag = setboundary(node,elem,'Dirichlet');
bdinfo = genMeshbdinfo(elem,bdFlag);
uclass = cell(7,9);
uclass{1,1} = @(x,y) (1-x.^2-y.^2)/4; % exact solution
bdinfo.function = cell(1,4);
bdinfo.function{1} = uclass{1,1};


%% FEM main routine
elemtype = size(mesh.T,1);

% generate trial and test basis
switch elemtype
    case 3  % triangular mesh
        basis.trial = genTribasis2(trbasis_type,trbasis_order);
        basis.test = genTribasis2(tsbasis_type,tsbasis_order);
    case 4  % quadrilateral mesh
        basis.trial = genQuadbasis2(trbasis_type,trbasis_order);
        basis.test = genQuadbasis2(tsbasis_type,tsbasis_order);
end

% generate Pb and Tb information matrix
trPT = genPbTb(mesh,trbasis_type,trbasis_order,femtype);
tsPT = genPbTb(mesh,tsbasis_type,tsbasis_order,femtype);
basis.trial.Pb = trPT.Pb; basis.trial.Tb = trPT.Tb;
basis.test.Pb = tsPT.Pb; basis.test.Tb = tsPT.Tb;
basis.trial.Nb = size(trPT.Pb,2); basis.trial.Nlb = size(trPT.Tb,1);
basis.test.Nb = size(tsPT.Pb,2); basis.test.Nlb = size(tsPT.Tb,1);

% assemble stiffness matrix and RHS matrix
const = @(x,y) 1;
A = assembTristiff(mesh,basis,const,1,0,1,0,Gauss_order,femtype,1);
A = A + assembTristiff(mesh,basis,const,0,1,0,1,Gauss_order,femtype,1);
B = assembTristiff(mesh,basis,const,0,0,0,0,Gauss_order,femtype,1);

% generate boundary information
fembdinfo.trial = genFEMbdinfo(bdinfo,trbasis_type,trbasis_order,trPT,elemtype,femtype);
fembdinfo.test = genFEMbdinfo(bdinfo,tsbasis_type,tsbasis_order,tsPT,elemtype,femtype);
fembdinfo.function = bdinfo.function;

% apply Dirichlet zero boundary
Ndof = size(A,1);
Dnode = fembdinfo.test.nodes(2,fembdinfo.test.nodes(1,:)==1);
fixedNode = Dnode;
if ~isempty(Dnode)
    gD = fembdinfo.function{1};
    bdidx = zeros(Ndof,1);
    bdidx(Dnode) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    A = T*A + Tbd;
    B = T*B + Tbd;
end

% solve the eigenvalue problem
[U,V] = eigs(A,B,k,'smallestabs','Tolerance',1e-8);
l1 = length(find(abs(diag(V)-1)<=1e-3));
U = U(:,l1+1:end);
V = diag(V(l1+1:end,l1+1:end));
lnum = length(V);
toc;

% compute spectral error
Efem = zeros(1,lnum); uh = zeros(size(U,1),1);
for i = 1:lnum
    uh = uh + U(:,i)/V(i) * semi_error(mesh,basis,U(:,i),0,0,Gauss_order);
    Efem(i) = max( abs( uclass{1,1}(mesh.P(1,:),mesh.P(2,:))' - uh ) );
end

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
set(gcf,'unit','normalized','position',[0.2, 0.2, 0.8, 0.6]);
semilogy(1:lnum,Efem,'-*',... %1:lnum,EI,'-o',1:lnum,EIfem,'-.x',1:lnum,R,'--+',
         'linewidth',1); hold on;
str1 = ['Lagrange ',num2str(trbasis_order),' order element ',soblev,' error'];
str2 = strcat(['$||u-u_h||_{',soblev,'}$']);
% str3 = strcat(['$||u-I_h(u)||_{',soblev,'}$']);
% str4 = strcat(['$||I_h(u)-u_h||_{',soblev,'}$']);
str5 = ['$y=O(h^',num2str(ord),')$'];
title(str1);
legend(str2,'location','southwest','interpreter','latex');
xlabel('$Num\ of\ eigenvalue$','Interpreter','latex');
ylabel('Error','Interpreter','latex');
set(gca,'fontsize',18);

% for i = 1:20
%     figure(1)
%     subplot(4,5,i)
%     pdesurf(p,t,U(:,i));
%     view(0,90)
%     colormap parula;
%     colorbar
%     lk = strcat(['$\lambda','(',num2str(i),') = ',num2str(V(i,i)),'$']);
%     title(lk,'Interpreter','latex');
% end

% for i = 1:20
%     figure(2)
%     subplot(4,5,i)
%     pdesurf(p,t,u(:,i));
%     view(0,90)
%     colormap parula;
%     colorbar
%     lk = strcat(['$\lambda','(',num2str(i),') = ',num2str(l(i)),'$']);
%     title(lk,'Interpreter','latex');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subroutine to compute int_{\Omega} d^{s+t}uh/dx^{s}dy^{t}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = semi_error(mesh,basis,uh,s,t,Gauss_order)

trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;
[intpt,weight] = Gauss_int_tri_ref2(Gauss_order);

pn1 = mesh.P(:,mesh.T(1,:));
v21 = mesh.P(:,mesh.T(2,:))-pn1;
v31 = mesh.P(:,mesh.T(3,:))-pn1;
detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);

s = num2str(s);  t = num2str(t); 
switch strcat(s,t)
    case '00'  % L^2 error
        Uhat = zeros(length(weight),trnlb);
        for k = 1:trnlb
            Uhat(:,k) = trf{k,1,1}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uhat * uh(trtb);
        % Rmk: in the loop above it should be tr{k,1,1}(intpt(:,1),intpt(:,2)),
        % not tr{k,1,1}(xn(:,1),xn(:,2)), since dphi is defined on refference element!

    case '10' % ux error
        J22 = v31(2,:)./detJ;
        nJ21 = -v21(2,:)./detJ;
        Uxhat = zeros(length(weight),trnlb); 
        Uyhat = Uxhat;
        for k = 1:trnlb
            Uxhat(:,k) = trf{k,2,2}(intpt(:,1),intpt(:,2));
            Uyhat(:,k) = trf{k,1,2}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uxhat * ( repmat(J22,trnlb,1).*uh(trtb) ) + ...
               Uyhat * ( repmat(nJ21,trnlb,1).*uh(trtb) ) ;

    case '01' % uy error
        nJ12 = -v31(1,:)./detJ;
        J11 = v21(1,:)./detJ;
        Uxhat = zeros(length(weight),trnlb);
        Uyhat = Uxhat;
        for k = 1:trnlb
            Uxhat(:,k) = trf{k,2,2}(intpt(:,1),intpt(:,2));
            Uyhat(:,k) = trf{k,1,2}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uxhat * ( repmat(nJ12,trnlb,1).*uh(trtb) ) + ...
               Uyhat * ( repmat(J11,trnlb,1).*uh(trtb) ) ;

    case '20' % uxx error
        J22 = v31(2,:)./detJ;
        nJ21 = -v21(2,:)./detJ;
        cn1 = J22.^2;  cn2 = 2*J22.*nJ21;  cn3 = nJ21.^2;
        Uxxhat = zeros(length(weight),trnlb);
        Uxyhat = Uxxhat; 
        Uyyhat = Uxxhat; 
        for k = 1:trnlb
            Uxxhat(:,k) = trf{k,3,3}(intpt(:,1),intpt(:,2));
            Uxyhat(:,k) = trf{k,2,3}(intpt(:,1),intpt(:,2));
            Uyyhat(:,k) = trf{k,1,3}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uxxhat * ( repmat(cn1,trnlb,1).*uh(trtb) ) + ...
               Uxyhat * ( repmat(cn2,trnlb,1).*uh(trtb) ) + ...
               Uyyhat * ( repmat(cn3,trnlb,1).*uh(trtb) );

    case '11' % uxy error
        J22 = v31(2,:)./detJ;
        nJ21 = -v21(2,:)./detJ;
        nJ12 = -v31(1,:)./detJ;
        J11 = v21(1,:)./detJ;
        cn1 = nJ12.*J22;  cn2 = nJ12.*nJ21 + J11.*J22;  cn3 = J11.*nJ21;
        Uxxhat = zeros(length(weight),trnlb);
        Uxyhat = Uxxhat; 
        Uyyhat = Uxxhat; 
        for k = 1:trnlb
            Uxxhat(:,k) = trf{k,3,3}(intpt(:,1),intpt(:,2));
            Uxyhat(:,k) = trf{k,2,3}(intpt(:,1),intpt(:,2));
            Uyyhat(:,k) = trf{k,1,3}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uxxhat * ( repmat(cn1,trnlb,1).*uh(trtb) ) + ...
               Uxyhat * ( repmat(cn2,trnlb,1).*uh(trtb) ) + ...
               Uyyhat * ( repmat(cn3,trnlb,1).*uh(trtb) );

    case '02' % uyy error
        nJ12 = -v31(1,:)./detJ;
        J11 = v21(1,:)./detJ;
        cn1 = nJ12.^2;  cn2 = 2*nJ12.*J11;  cn3 = J11.^2;
        Uxxhat = zeros(length(weight),trnlb);
        Uxyhat = Uxxhat; 
        Uyyhat = Uxxhat; 
        for k = 1:trnlb
            Uxxhat(:,k) = trf{k,3,3}(intpt(:,1),intpt(:,2));
            Uxyhat(:,k) = trf{k,2,3}(intpt(:,1),intpt(:,2));
            Uyyhat(:,k) = trf{k,1,3}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uxxhat * ( repmat(cn1,trnlb,1).*uh(trtb) ) + ...
               Uxyhat * ( repmat(cn2,trnlb,1).*uh(trtb) ) + ...
               Uyyhat * ( repmat(cn3,trnlb,1).*uh(trtb) );
end

% area(ref_angle) = 1/2
y = sum( detJ/2.*( weight * (ipt2) ) );  
%y = abs(sqrt(y));

end




