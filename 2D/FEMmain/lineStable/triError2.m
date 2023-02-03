function y = triError2(mesh,basis,uclass,uh,errtype,Gauss_order,femtype)
% get error between exact solution u and FEM solution uh in triangluar mesh
% uclass : u{i,j}(x,y) = d^{i-1}u/dx^{j-1}dy^{i-j}(x,y) 
% uh : FEM solution, i.e. the coefficient of each global basis function
% errtype: 'linf' for L^infty norm error; 'l2' for L^2 norm error
%       'h1' for H^1 norm error; 'h1semi' for H^1 semi-norm error, etc.

switch errtype
    case 'linf' % choose Gauss quadrature points in each element to estimate max error
        [intpt, ~] = Gauss_int_tri_ref2(Gauss_order);
        v21 = mesh.P(:,mesh.T(2,:))-mesh.P(:,mesh.T(1,:));
        v31 = mesh.P(:,mesh.T(3,:))-mesh.P(:,mesh.T(1,:));
        J = [reshape(v21,[],1),reshape(v31,[],1)];
        xn = J*intpt'+reshape(mesh.P(:,mesh.T(1,:)),[],1);
        xn = xn';
        pxy = reshape(xn,[],mesh.N);
        px = reshape(pxy(1:length(intpt(:,1)),:),[],1);
        py = reshape(pxy(length(intpt(:,1))+1:end,:),[],1);
        ipt1 = reshape(uclass{1,1}(px,py),[],mesh.N);  % exact value on Gauss quadrature
        clear pxy xn px py
        Uhat = zeros(length(intpt(:,1)),basis.trial.Nlb);
        for k = 1:basis.trial.Nlb
            Uhat(:,k) = basis.trial.reffun{k}(intpt(:,1),intpt(:,2));
        end
        ipt2 = Uhat * uh(basis.trial.Tb);
        y = max(max(abs(ipt1-ipt2)));
    case 'l2'
        y = semi_error(mesh,basis,uclass{1,1},uh,0,0,Gauss_order);
    case 'h1semi'
        y1 = semi_error(mesh,basis,uclass{2,2},uh,1,0,Gauss_order); % ux
        y2 = semi_error(mesh,basis,uclass{2,1},uh,0,1,Gauss_order); % uy
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h1'
        y1 = triError2(mesh,basis,uclass,uh,'l2',Gauss_order,femtype);
        y2 = triError2(mesh,basis,uclass,uh,'h1semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h2semi'
        y1 = semi_error(mesh,basis,uclass{3,3},uh,2,0,Gauss_order); % uxx
        y2 = semi_error(mesh,basis,uclass{3,2},uh,1,1,Gauss_order); % uxy
        y3 = semi_error(mesh,basis,uclass{3,1},uh,0,2,Gauss_order); % uyy
        y = abs( sqrt( y1^2 + y2^2 + y3^2 ) );
    case 'h2'
        y1 = triError2(mesh,basis,uclass,uh,'h1',Gauss_order,femtype);
        y2 = triError2(mesh,basis,uclass,uh,'h2semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h3semi'
        y = zeros(1,4);
        for i = 1:4
            y(i) = semi_error(mesh,basis,uclass{4,i},uh,4-i,i-1,Gauss_order);
        end
        y = norm(y);
    case 'h3'
        y1 = triError2(mesh,basis,uclass,uh,'h2',Gauss_order,femtype);
        y2 = triError2(mesh,basis,uclass,uh,'h3semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h4semi'
        y = zeros(1,5);
        for i = 1:5
            y(i) = semi_error(mesh,basis,uclass{5,i},uh,5-i,i-1,Gauss_order);
        end
        y = norm(y);
    case 'h4'
        y1 = triError2(mesh,basis,uclass,uh,'h3',Gauss_order,femtype);
        y2 = triError2(mesh,basis,uclass,uh,'h4semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h5semi'
        y = zeros(1,6);
        for i = 1:6
            y(i) = semi_error(mesh,basis,uclass{6,i},uh,6-i,i-1,Gauss_order);
        end
        y = norm(y);
    case 'h5'
        y1 = triError2(mesh,basis,uclass,uh,'h4',Gauss_order,femtype);
        y2 = triError2(mesh,basis,uclass,uh,'h5semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h6semi'
        y = zeros(1,7);
        for i = 1:7
            y(i) = semi_error(mesh,basis,uclass{7,i},uh,7-i,i-1,Gauss_order);
        end
        y = norm(y);
    case 'h6'
        y1 = triError2(mesh,basis,uclass,uh,'h5',Gauss_order,femtype);
        y2 = triError2(mesh,basis,uclass,uh,'h6semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
end

% if strcmp(femtype.type,'ipdg') == 1
%     switch type
%         case 'ipdg_h1'
%             y = semi_error(mesh,basis,uclass{1},uh,0,Gauss_order);
%             
%     end
% end

end



function y = semi_error(mesh,basis,du,uh,s,t,Gauss_order)
% input: d^(s+t)u/dx^sdy^t and FEM solution uh, 
% output: |d^(s+t)(u - uh)/dx^sdy^t|_{L^2}

trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;
[intpt,weight] = Gauss_int_tri_ref2(Gauss_order);

pn1 = mesh.P(:,mesh.T(1,:));
v21 = mesh.P(:,mesh.T(2,:))-pn1;
v31 = mesh.P(:,mesh.T(3,:))-pn1;
detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);
J = [reshape(v21,[],1),reshape(v31,[],1)];
xn = J*intpt'+reshape(pn1,[],1);
xn = xn';
pxy = reshape(xn,[],mesh.N);
px = reshape(pxy(1:length(weight),:),[],1);
py = reshape(pxy(length(weight)+1:end,:),[],1);

ipt1 = reshape(du(px,py),[],mesh.N);  % exact value on Gauss quadrature

clear pxy xn px py
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
y = sum( detJ/2.*( weight * (ipt1 - ipt2).^2 ) );  
y = abs(sqrt(y));

end


