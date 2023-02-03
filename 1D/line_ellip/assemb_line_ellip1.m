function [A,F] = assemb_line_ellip1(mesh,basis,coeff,f,bdinfo,Gauss_nodes,femtype)
% assemble stiffness matrix and load vector£¬then deal with boundary condition
cfc = coeff.c; 
cfd = coeff.d;
cfk = coeff.k;
trnb = basis.trial.Nb; trnlb = basis.trial.Nlb;
tsnb = basis.test.Nb; tsnlb = basis.test.Nlb;
trf = basis.trial.reffun;
tsf = basis.test.reffun;
trtb = basis.trial.Tb;
tstb = basis.test.Tb;
trpb = basis.trial.Pb;
tspb = basis.test.Pb;
[intpt,weight] = Gauss_int_ref1(Gauss_nodes); 
eps = 1e-14;
%S = zeros(tsnlb,trnlb);
tr = cell(2,trnlb);
ts = cell(2,tsnlb);
for i=1:trnlb
    for j=1:tsnlb
        tr{2,i} = str2func(func2str(trf{2,i}));
        tr{1,i} = str2func(func2str(trf{1,i}));
        ts{2,j} = str2func(func2str(tsf{2,j}));
        ts{1,j} = str2func(func2str(tsf{1,j}));
    end
end

% assemble stiffness matrix A
A = sparse(tsnb,trnb);
switch femtype.type
    case 'pgm' % PG method
        % row1 = []; col1 = []; s1 = [];  % method 3
        % parallel computing
        for n=1:mesh.N
            h = mesh.P(:,mesh.T(2,n))-mesh.P(:,mesh.T(1,n));
            w = (mesh.P(:,mesh.T(1,n))+mesh.P(:,mesh.T(2,n)))/2;
            row2 = zeros(1,trnlb*tsnlb); col2 = row2; s2 = row2; 
            count = 1;
            for i=1:trnlb
                for j=1:tsnlb
                    % assemble element stiffness matrix
                    % each derivative of basis function contributes a factor of 2/h
                    ipt = cfc(h/2*intpt+w).*tr{2,i}(intpt).*ts{2,j}(intpt)*(2/h)^2 ...
                         + cfd(h/2*intpt+w).*tr{2,i}(intpt).*ts{1,j}(intpt).*(2/h) ...
                         + cfk(h/2*intpt+w).*tr{1,i}(intpt).*ts{1,j}(intpt);
                    row2(count) = tstb(j,n);
                    col2(count) = trtb(i,n);  
                    s2(count) = h/2 * weight * ipt';
                    count = count + 1;
                    %S(j,i) = h/2 * weight * ipt';
                end
            end
            % add element stiffness matrix to global stiffness matrix
            % method 1:
            %      A(tstb(:,n),trtb(:,n)) = ...
            %      A(tstb(:,n),trtb(:,n)) + S;
            % method 2:
                 A = A + sparse(row2,col2,s2,tsnb,trnb);
            % method 3:
%                  [row1,col1,s1] = concat1(row1,col1,s1,row2,col2,s2,...
%                                          trnlb,tsnlb,ordtr,ordts);
        end
%         A = A + sparse(row1,col1,s1,tsnb,trnb); % method 3
        
    case 'ipdg' % IPDG method
        rr1 = []; cc1 = []; ss1 = [];
        rr2 = []; cc2 = []; ss2 = [];
        rr3 = []; cc3 = []; ss3 = [];
        for n=1:mesh.N
            h = mesh.P(:,mesh.T(2,n))-mesh.P(:,mesh.T(1,n));
            w = (mesh.P(:,mesh.T(1,n))+mesh.P(:,mesh.T(2,n)))/2;
            row1 = zeros(1,trnlb*tsnlb); col1 = row1; s1 = row1; 
            [row2,col2,s2] = find(zeros(tsnlb,trnlb));
            [row3,col3,s3] = find(zeros(tsnlb,trnlb));
            count = 1;
            for i=1:trnlb
                for j=1:tsnlb
                    % assemble element stiffness matrix on the main diagonal
                    % each derivative of basis function contributes a factor of 2/h 
                    ipt = cfc(h/2*intpt+w).*tr{2,i}(intpt).*ts{2,j}(intpt)*(2/h)^2 ...
                         + cfd(h/2*intpt+w).*tr{2,i}(intpt).*ts{1,j}(intpt).*(2/h) ...
                         + cfk(h/2*intpt+w).*tr{1,i}(intpt).*ts{1,j}(intpt);
                    row1(count) = tstb(j,n);
                    col1(count) = trtb(i,n);
                    s1(count) = h/2 * weight * ipt';
 
                    % assemble element stiffness matrix, main and sub diagonal 
                    if n==1
                        w2 = (mesh.P(:,mesh.T(1,n+1))+mesh.P(:,mesh.T(2,n+1)))/2;
                        % trial and test function both in 1st element
                        ptr = mesh.P(:,mesh.T(2,n))-eps;
                        hr = mesh.P(:,mesh.T(2,n+1))-mesh.P(:,mesh.T(1,n+1));
                        he = max(h,hr);
                        cdu0 = 1/2*cfc(ptr)*tr{2,i}(2/h*(ptr-w))*(2/h);
                        cdv0 = 1/2*cfc(ptr)*ts{2,j}(2/h*(ptr-w))*(2/h);
                        cdu1 = cfc(ptr)*tr{2,i}(2/h*(ptr-w))*(2/h);
                        cdv1 = cfc(ptr)*ts{2,j}(2/h*(ptr-w))*(2/h);
                        u1 = tr{1,i}(2/h*(ptr-w));
                        v1 = ts{1,j}(2/h*(ptr-w));
                        s1(count) = s1(count) - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                            + femtype.ipdg.gamma0(mesh.T(2,n))*he^(-femtype.ipdg.alpha(mesh.T(2,n)))*u1*v1...
                            + femtype.ipdg.gamma1(mesh.T(2,n))*he^(-femtype.ipdg.beta(mesh.T(2,n)))*cdu1*cdv1;
                        
                        % trial in 1st element and test in 2nd element
                        ptr1 = mesh.P(:,mesh.T(2,n))-eps;
                        ptr2 = mesh.P(:,mesh.T(2,n))+eps;
                        hr = mesh.P(:,mesh.T(2,n+1))-mesh.P(:,mesh.T(1,n+1));
                        he = max(h,hr);
                        cdu0 = 1/2*cfc(ptr1)*tr{2,i}(2/h*(ptr1-w))*(2/h);
                        cdv0 = 1/2*cfc(ptr2)*ts{2,j}(2/hr*(ptr2-w2))*(2/hr);
                        cdu1 = cfc(ptr1)*tr{2,i}(2/h*(ptr1-w))*(2/h);
                        cdv1 = -cfc(ptr2)*ts{2,j}(2/hr*(ptr2-w2))*(2/hr);  
                        u1 = tr{1,i}(2/h*(ptr1-w));
                        v1 = -ts{1,j}(2/hr*(ptr2-w2));
                        row2(count) = tstb(j,n)+tsnlb;
                        col2(count) = trtb(i,n);
                        s2(count) = - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                             + femtype.ipdg.gamma0(mesh.T(2,n))*he^(-femtype.ipdg.alpha(mesh.T(2,n)))*u1*v1...
                             + femtype.ipdg.gamma1(mesh.T(2,n))*he^(-femtype.ipdg.beta(mesh.T(2,n)))*cdu1*cdv1;
                         
                    elseif n==mesh.N
                        w1 = (mesh.P(:,mesh.T(1,n-1))+mesh.P(:,mesh.T(2,n-1)))/2;
                        % trial and test function both in Nth element
                        ptl = mesh.P(:,mesh.T(1,n))+eps;
                        hl = mesh.P(:,mesh.T(2,n-1))-mesh.P(:,mesh.T(1,n-1));
                        he = max(h,hl);
                        cdu0 = 1/2*cfc(ptl)*tr{2,i}(2/h*(ptl-w))*(2/h);
                        cdv0 = 1/2*cfc(ptl)*ts{2,j}(2/h*(ptl-w))*(2/h);
                        cdu1 = -cfc(ptl)*tr{2,i}(2/h*(ptl-w))*(2/h);
                        cdv1 = -cfc(ptl)*ts{2,j}(2/h*(ptl-w))*(2/h);
                        u1 = -tr{1,i}(2/h*(ptl-w));
                        v1 = -ts{1,j}(2/h*(ptl-w));
                        s1(count) = s1(count) - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                            + femtype.ipdg.gamma0(mesh.T(1,n))*he^(-femtype.ipdg.alpha(mesh.T(1,n)))*u1*v1...
                            + femtype.ipdg.gamma1(mesh.T(1,n))*he^(-femtype.ipdg.beta(mesh.T(1,n)))*cdu1*cdv1;
                        
                        % trial in Nth element and test in (N-1)th element
                        ptl1 = mesh.P(:,mesh.T(1,n))-eps;
                        ptl2 = mesh.P(:,mesh.T(1,n))+eps;
                        hl = mesh.P(:,mesh.T(2,n-1))-mesh.P(:,mesh.T(1,n-1));
                        he = max(h,hl);
                        cdu0 = 1/2*cfc(ptl2)*tr{2,i}(2/h*(ptl2-w))*(2/h);
                        cdv0 = 1/2*cfc(ptl1)*ts{2,j}(2/hl*(ptl1-w1))*(2/hl);
                        cdu1 = -cfc(ptl2)*tr{2,i}(2/h*(ptl2-w))*(2/h);
                        cdv1 = cfc(ptl1)*ts{2,j}(2/hl*(ptl1-w1))*(2/hl);  
                        u1 = -tr{1,i}(2/h*(ptl2-w));
                        v1 = ts{1,j}(2/hl*(ptl1-w1));
                        row3(count) = tstb(j,n)-tsnlb;
                        col3(count) = trtb(i,n);
                        s3(count) = - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                             + femtype.ipdg.gamma0(mesh.T(1,n))*he^(-femtype.ipdg.alpha(mesh.T(1,n)))*u1*v1...
                             + femtype.ipdg.gamma1(mesh.T(1,n))*he^(-femtype.ipdg.beta(mesh.T(1,n)))*cdu1*cdv1;
                    else
                        w1 = (mesh.P(:,mesh.T(1,n-1))+mesh.P(:,mesh.T(2,n-1)))/2;
                        w2 = (mesh.P(:,mesh.T(1,n+1))+mesh.P(:,mesh.T(2,n+1)))/2;
                        % right node
                        % trial and test function both in nth element
                        ptr = mesh.P(:,mesh.T(2,n))-eps;
                        hr = mesh.P(:,mesh.T(2,n+1))-mesh.P(:,mesh.T(1,n+1));
                        he = max(h,hr);
                        cdu0 = 1/2*cfc(ptr)*tr{2,i}(2/h*(ptr-w))*(2/h);
                        cdv0 = 1/2*cfc(ptr)*ts{2,j}(2/h*(ptr-w))*(2/h);
                        cdu1 = cfc(ptr)*tr{2,i}(2/h*(ptr-w))*(2/h);
                        cdv1 = cfc(ptr)*ts{2,j}(2/h*(ptr-w))*(2/h);
                        u1 = tr{1,i}(2/h*(ptr-w));
                        v1 = ts{1,j}(2/h*(ptr-w));
                        s1(count) = s1(count) - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                            + femtype.ipdg.gamma0(mesh.T(2,n))*he^(-femtype.ipdg.alpha(mesh.T(2,n)))*u1*v1...
                            + femtype.ipdg.gamma1(mesh.T(2,n))*he^(-femtype.ipdg.beta(mesh.T(2,n)))*cdu1*cdv1;
                        
                        % trial in nth element and test in (n+1)th element
                        ptr1 = mesh.P(:,mesh.T(2,n))-1e-14;
                        ptr2 = mesh.P(:,mesh.T(2,n))+1e-14;
                        hr = mesh.P(:,mesh.T(2,n+1))-mesh.P(:,mesh.T(1,n+1));
                        he = max(h,hr);
                        cdu0 = 1/2*cfc(ptr1)*tr{2,i}(2/h*(ptr1-w))*(2/h);
                        cdv0 = 1/2*cfc(ptr2)*ts{2,j}(2/hr*(ptr2-w2))*(2/hr);
                        cdu1 = cfc(ptr1)*tr{2,i}(2/h*(ptr1-w))*(2/h);
                        cdv1 = -cfc(ptr2)*ts{2,j}(2/hr*(ptr2-w2))*(2/hr);  
                        u1 = tr{1,i}(2/h*(ptr1-w));
                        v1 = -ts{1,j}(2/hr*(ptr2-w2));
                        row2(count) = tstb(j,n)+tsnlb;
                        col2(count) = trtb(i,n);
                        s2(count) = - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                             + femtype.ipdg.gamma0(mesh.T(2,n))*he^(-femtype.ipdg.alpha(mesh.T(2,n)))*u1*v1...
                             + femtype.ipdg.gamma1(mesh.T(2,n))*he^(-femtype.ipdg.beta(mesh.T(2,n)))*cdu1*cdv1;
                        
                        % left node
                        % trial and test function both in nth element
                        ptl = mesh.P(:,mesh.T(1,n))+eps;
                        hl = mesh.P(:,mesh.T(2,n-1))-mesh.P(:,mesh.T(1,n-1));
                        he = max(h,hl);
                        cdu0 = 1/2*cfc(ptl)*tr{2,i}(2/h*(ptl-w))*(2/h);
                        cdv0 = 1/2*cfc(ptl)*ts{2,j}(2/h*(ptl-w))*(2/h);
                        cdu1 = -cfc(ptl)*tr{2,i}(2/h*(ptl-w))*(2/h);
                        cdv1 = -cfc(ptl)*ts{2,j}(2/h*(ptl-w))*(2/h);
                        u1 = -tr{1,i}(2/h*(ptl-w));
                        v1 = -ts{1,j}(2/h*(ptl-w));
                        s1(count) = s1(count) - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                            + femtype.ipdg.gamma0(mesh.T(1,n))*he^(-femtype.ipdg.alpha(mesh.T(1,n)))*u1*v1...
                            + femtype.ipdg.gamma1(mesh.T(1,n))*he^(-femtype.ipdg.beta(mesh.T(1,n)))*cdu1*cdv1;
                        
                        % trial in nth element and test in (n-1)th element
                        ptl1 = mesh.P(:,mesh.T(1,n))-eps;
                        ptl2 = mesh.P(:,mesh.T(1,n))+eps;
                        hl = mesh.P(:,mesh.T(2,n-1))-mesh.P(:,mesh.T(1,n-1));
                        he = max(h,hl);
                        cdu0 = 1/2*cfc(ptl2)*tr{2,i}(2/h*(ptl2-w))*(2/h);
                        cdv0 = 1/2*cfc(ptl1)*ts{2,j}(2/hl*(ptl1-w1))*(2/hl);
                        cdu1 = -cfc(ptl2)*tr{2,i}(2/h*(ptl2-w))*(2/h);
                        cdv1 = cfc(ptl1)*ts{2,j}(2/hl*(ptl1-w1))*(2/hl);  
                        u1 = -tr{1,i}(2/h*(ptl2-w));
                        v1 = ts{1,j}(2/hl*(ptl1-w1));
                        row3(count) = tstb(j,n)-tsnlb;
                        col3(count) = trtb(i,n);
                        s3(count) = - cdu0 * v1 - femtype.ipdg.theta * cdv0 * u1 ...
                             + femtype.ipdg.gamma0(mesh.T(1,n))*he^(-femtype.ipdg.alpha(mesh.T(1,n)))*u1*v1...
                             + femtype.ipdg.gamma1(mesh.T(1,n))*he^(-femtype.ipdg.beta(mesh.T(1,n)))*cdu1*cdv1;
                        
                    end
                    count = count + 1;
                end
            end
            % add element stiffness matrix to global stiffness matrix
            % method 1:
%             A = A + sparse(row1,col1,s1,tsnb,trnb) + sparse(row2,col2,s2,tsnb,trnb) ...
%                 + sparse(row3,col3,s3,tsnb,trnb);
            % method 2:
            rr1 = [rr1, row1]; cc1 = [cc1, col1]; ss1 = [ss1, s1];
            rr2 = [rr2, row2]; cc2 = [cc2, col2]; ss2 = [ss2, s2];
            rr3 = [rr3, row3]; cc3 = [cc3, col3]; ss3 = [ss3, s3];
        end
        A = A + sparse(rr1,cc1,ss1,tsnb,trnb) + sparse(rr2,cc2,ss2,tsnb,trnb) ...
                + sparse(rr3,cc3,ss3,tsnb,trnb);
end

% assemble load vector F
F = zeros(tsnb,1);
d = zeros(tsnlb,1);
for n=1:mesh.N 
    h = mesh.P(:,mesh.T(2,n))-mesh.P(:,mesh.T(1,n));
    w = (mesh.P(:,mesh.T(1,n))+mesh.P(:,mesh.T(2,n)))/2;
        for j=1:tsnlb
            ipt = f(h/2*intpt+w).*ts{1,j}(intpt);
            d(j,1) = h/2 * weight * ipt';
        end
    F(tstb(:,n),1) = ...
    F(tstb(:,n),1) + d;
end

% add boundary condition, b.c. 
for k=1:size(bdinfo.nodes,2)
    mesh.i = bdinfo.nodes(2,k);
    tmp = find(tspb==mesh.P(:,mesh.i));
    i = tmp(1);
    tmp = find(trpb==mesh.P(:,mesh.i));
    j = tmp(1);
    if bdinfo.nodes(1,k)==1   % Dirichlet boundary condition
        A(i,:)=0;
        A(i,j)=1;
        F(i)=bdinfo.value(1,k);
    elseif bdinfo.nodes(1,k)==2   % Neumann boundary condition
        if i==1   % Neumann b.c. on left node a
           F(i)=F(i)-bdinfo.value(1,k)*cfc(mesh.left);
        else      % Neumann b.c. on right node b
           F(i)=F(i)+bdinfo.value(1,k)*cfc(mesh.right); 
        end
    elseif bdinfo.nodes(1,k)==3   % Robin boundary condition
        if i==1   % Robin b.c. on left node a
           A(i,j)=A(i,j)-bdinfo.value(1,k)*cfc(mesh.left);
           F(i)=F(i)-bdinfo.value(2,k)*cfc(mesh.left);
        else      % Robin b.c. on right node b
           A(i,j)=A(i,j)+bdinfo.value(1,k)*cfc(mesh.right);
           F(i)=F(i)+bdinfo.value(2,k)*cfc(mesh.right); 
        end
    else
        error('bdinfo.nodes(1,:) must be 1, 2 or 3');
    end
end

pure = strcmp(func2str(cfk),'@(x)0');
if bdinfo.nodes(1,1)==2 && bdinfo.nodes(1,2)==2 && pure==1 
% pure Nuemann b.c., add extra condition: int_a^b u_j phi_j = 0
    B = sparse(tsnb+1,trnb);
    B(1:end-1,:) = A;
    d = zeros(1,trnlb);
    for n=1:mesh.N 
    h = mesh.P(:,mesh.T(2,n))-mesh.P(:,mesh.T(1,n));
        for j=1:trnlb
            d(1,j) = h/2 * weight * tr{1,j}(intpt)';
        end
    B(end,trtb(:,n)) = B(end,trtb(:,n)) + d;
    end
    A = B;
    G = zeros(trnb+1,1);
    G(1:end-1,1) = F;
    F = G;
end


end

