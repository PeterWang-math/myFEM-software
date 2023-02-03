function A = assembstiff(mesh,basis,coeff,Gauss_order,femtype)
% function A = assembstiff(mesh,basis,c,s,t,p,q,Gauss_order,femtype)
% assemble stiffness matrix A for 2d linear ellipse case
% assemble int_{En} c x d^{s+t}phi/dx^{s}dy^{t} x d^{p+q}psai/dx^{p}dy^{q}

cfc = coeff.c; cnum = numel(cfc);
cfd = coeff.d;
cfk = coeff.k;
trnb = basis.trial.Nb; trnlb = basis.trial.Nlb;
tsnb = basis.test.Nb; tsnlb = basis.test.Nlb;
trf = basis.trial.reffun;
tsf = basis.test.reffun;

switch size(mesh.T,1)
    case 3
        [intpt,weight] = Gauss_int_tri_ref2(Gauss_order); 
        eleArea = 1/2;
    case 4
        [intpt,weight] = Gauss_int_quad_ref2(Gauss_order); 
        eleArea = 1;
end


% assemble stiffness matrix A
% 去掉关于单元N的循环
A = sparse(tsnb,trnb);
switch femtype.type
    case 'pgm' % PG method
        switch eleArea
            case 1/2  % triangle mesh
                pn1 = mesh.P(:,mesh.T(1,:));
                v21 = mesh.P(:,mesh.T(2,:))-pn1;
                v31 = mesh.P(:,mesh.T(3,:))-pn1;
                detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);
                J = [reshape(v21,[],1),reshape(v31,[],1)];
            case 1  % quadrilateral mesh
                pn1 = mesh.P(:,mesh.T(1,:));
                pn2 = mesh.P(:,mesh.T(2,:));
                pn3 = mesh.P(:,mesh.T(3,:));
                pn4 = mesh.P(:,mesh.T(4,:));
                v21 = pn2-pn1;
                v41 = pn4-pn1;
                detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);
                J = [abs(mesh.P(1,mesh.T(4,n))-mesh.P(1,mesh.T(1,n))),0;
                    0,abs(mesh.P(2,mesh.T(2,n))-mesh.P(2,mesh.T(1,n)))];
        end

        xn = J*intpt'+reshape(pn1,[],1);
        xn = xn';
        pxy = reshape(xn,[],mesh.N);
        px = reshape(pxy(1:length(weight),:),[],1);
        py = reshape(pxy(length(weight)+1:end,:),[],1);
        clear pxy xn
        
        % assemble element stiffness matrix
        for i=1:tsnlb
            for j=1:trnlb
                % compute dphi/dx by chain rule
                %% ts for i && tr for j ! ! !  每行是test，每列是trial
                switch eleArea
                    case 1/2    % triangle mesh
                        dtrdx = trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                            trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                        dtrdy =-trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                            trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                        dtsdx = tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                            tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                        dtsdy =-tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                            tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                    case 1    % quadrilateral mesh
                        dtrdx = trf{j,1,2}(intpt(:,1),intpt(:,2))/J(1,1);
                        dtrdy = trf{j,2,2}(intpt(:,1),intpt(:,2))/J(2,2);
                        dtsdx = tsf{i,1,2}(intpt(:,1),intpt(:,2))/J(1,1);
                        dtsdy = tsf{i,2,2}(intpt(:,1),intpt(:,2))/J(2,2);
                end
                dtrdx = reshape(dtrdx,[],1);
                dtrdy = reshape(dtrdy,[],1);
                dtsdx = reshape(dtsdx,[],1);
                dtsdy = reshape(dtsdy,[],1);
                switch cnum
                    case 1  % c is a constant
                        ipt1 = cfc(px,py).*( dtrdx.*dtsdx + dtrdy.*dtsdy );
                    case 4  % c is a 2x2 cell
                        ipt1 = ( cfc{1,1}(px,py).*dtrdx + cfc{1,2}(px,py).*dtrdy ).*dtsdx ...
                            + ( cfc{2,1}(px,py).*dtrdx + cfc{2,2}(px,py).*dtrdy ).*dtsdy ;
                end
                % d is a 1x2 cell
                ipt2 = ( cfd{1}(px,py).*dtrdx + cfd{2}(px,py).*dtrdy ).* ...
                        repmat(tsf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                % k is a scalar
                ipt3 = cfk(px,py).*repmat(trf{j,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1).*  ...
                        repmat(tsf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                % |hat{K}| = 1/2 for triangle and 1 for quadrilateral, must multiply it
                Aij = eleArea * detJ.*( weight * reshape(ipt1 + ipt2 + ipt3,[],mesh.N) );
                
                % add element stiffness matrix to global stiffness matrix
                A = A + sparse(basis.test.Tb(i,:),basis.trial.Tb(j,:),Aij,tsnb,trnb);
            end
        end
end
%clear dtrdx dtrdy dtsdx dtsdy ipt1 ipt2 ipt3 Aij

end
