function A = assembTristiff(mesh,basis,C,s,t,p,q,Gauss_order,femtype)
% assemble stiffness matrix A for 2d linear ellipse case
% assemble int_{En} C x d^{s+t}phi/dx^{s}dy^{t} x d^{p+q}psai/dx^{p}dy^{q}, C is a constant or function
% issymm: the stiffness matrix A is symmetric or not
% issame: the basis.trial.reffun = basis.test.reffun or not
%         when issame = 0, C must be a function!

% to be removed if Lagrange quadratic element for constant C is completed  
if isnumeric(C) && basis.trial.order >= 2 && basis.test.order >= 2
    C = str2func(strcat(['@(x,y)',num2str(C)]));
end

issame = 1;
if isnumeric(C) && ~issame
    warning('basis.trial.reffun ~= basis.test.reffun, please reset c as a function');
    exit;
end

if basis.trial.type == basis.test.type && basis.trial.order == basis.test.order ...
   && s == p && t == q
    issymm = 1;
else
    issymm = 0;
end

trnb = basis.trial.Nb; trnlb = basis.trial.Nlb;
tsnb = basis.test.Nb; tsnlb = basis.test.Nlb;
trf = basis.trial.reffun;
tsf = basis.test.reffun;
trtb = basis.trial.Tb;
tstb = basis.test.Tb;

A = sparse(tsnb,trnb); % Initialize A = 0

if isreal(C) && (C==0) % C = 0 
    return
end
stpq = strcat([num2str(s),num2str(t),num2str(p),num2str(q)]);
[intpt,weight] = Gauss_int_tri_ref2(Gauss_order);
eleArea = 1/2;

% assemble stiffness matrix A
switch femtype.type
    case 'pgm' % PG method
        pn1 = mesh.P(:,mesh.T(1,:));
        v21 = mesh.P(:,mesh.T(2,:))-pn1;
        v31 = mesh.P(:,mesh.T(3,:))-pn1;
        detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);
        J = [reshape(v21,[],1),reshape(v31,[],1)];

        if ~isnumeric(C) || ~(basis.trial.type == basis.test.type ...
                           && basis.trial.order == basis.test.order)
        % C is a function, use Gauss quadrature to compute integral
            xn = J*intpt'+reshape(pn1,[],1);
            xn = xn';
            pxy = reshape(xn,[],mesh.N);
            px = reshape(pxy(1:length(weight),:),[],1);
            py = reshape(pxy(length(weight)+1:end,:),[],1);
            clear pxy xn

            % assemble element stiffness matrix
            % compute dphi/dx by chain rule
            %% ts for i && tr for j ! ! ! 
            for i = 1:tsnlb
                for j = 1:trnlb
                    if issymm && i>j  % for symmetric matrix only computes upper triangular parts
                        continue;
                    end
                    switch stpq
                        case '0000'
                            ipt = C(px,py).*repmat(trf{j,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1).*  ...
                                   repmat(tsf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );    
                        case '1000'
                            dtrdx = trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                                    trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                            dtrdx = reshape(dtrdx,[],1);
                            ipt = C(px,py).*dtrdx.*repmat(tsf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );   
                        case '0100'
                            dtrdy = -trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                                    trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                            dtrdy = reshape(dtrdy,[],1);
                            ipt = C(px,py).*dtrdy.*repmat(tsf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '0010'
                            dtsdx = tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                                    tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                            dtsdx = reshape(dtsdx,[],1);
                            ipt = C(px,py).*dtsdx.*repmat(trf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '0001'
                            dtsdy = -tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                                    tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                            dtsdy = reshape(dtsdy,[],1);
                            ipt = C(px,py).*dtsdy.*repmat(trf{i,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '1010'
                            dtrdx = trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                                    trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                            dtsdx = tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                                    tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                            dtrdx = reshape(dtrdx,[],1);
                            dtsdx = reshape(dtsdx,[],1);
                            ipt = C(px,py).*dtrdx.*dtsdx;
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '1001'
                            dtrdx = trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                                    trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                            dtsdy = -tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                                    tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                            dtrdx = reshape(dtrdx,[],1);
                            dtsdy = reshape(dtsdy,[],1);
                            ipt = C(px,py).*dtrdx.*dtsdy;
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '0110'
                            dtrdy = -trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                                    trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                            dtsdx = tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(2,:)./detJ) - ...
                                    tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(2,:)./detJ);
                            dtrdy = reshape(dtrdy,[],1);
                            dtsdx = reshape(dtsdx,[],1);
                            ipt = C(px,py).*dtrdy.*dtsdx;
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '0101'
                            dtrdy = -trf{j,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                                    trf{j,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                            dtsdy = -tsf{i,2,2}(intpt(:,1),intpt(:,2))*(v31(1,:)./detJ) + ...
                                    tsf{i,1,2}(intpt(:,1),intpt(:,2))*(v21(1,:)./detJ);
                            dtrdy = reshape(dtrdy,[],1);
                            dtsdy = reshape(dtsdy,[],1);
                            ipt = C(px,py).*dtrdy.*dtsdy;
                            Aij = eleArea * detJ.*( weight * reshape(ipt,[],mesh.N) );
                        case '2000'

                        case '0200'

                        case '0020'

                        case '0002'

                        case '1020'

                        case '1011'

                        case '1002'

                        case '0120'

                        case '0111'

                        case '0102'

                        case '2010'

                        case '2001'

                        case '1110'

                        case '1101'

                        case '0210'

                        case '0201'

                        case '2020'

                        case '2002'

                        case '0220'

                        case '0202'

                    end
                    
                    % add element stiffness matrix to global stiffness matrix
                    if issymm
                        if (j==i)
                            A = A + sparse(tstb(i,:),trtb(j,:),Aij,tsnb,trnb);
                        else
                            A = A + sparse([tstb(i,:), tstb(j,:)],...
                                    [trtb(j,:), trtb(i,:)],[Aij, Aij],tsnb,trnb);
                        end
                    else
                        A = A + sparse(tstb(i,:), trtb(j,:),Aij,tsnb,trnb);
                    end

                end
            end
        else  % C is a constant, use barycentric coordinates formula to compute integral
            switch basis.test.type
                case 'l'
                    switch basis.test.order
                        case 1  % Lagrange linear basis
                            J22 = v31(2,:)./detJ;
                            nJ21 = -v21(2,:)./detJ;
                            nJ12 = -v31(1,:)./detJ;
                            J11 = v21(1,:)./detJ;
                            c1 = J22;  c2 = nJ21;  ct1 = nJ12;  ct2 = J11;
                            switch stpq
                                case '0000'
                                    A0 = eleArea * detJ * C * 1/12;
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0,A0,A0],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(1,:),tstb(2,:)],...
                                                   [trtb(2,:),trtb(3,:),trtb(3,:)],...
                                                   [A0,A0,A0],tsnb,trnb);
                                    A = A + A';
                                case '1000'
                                    Ax01 = eleArea * detJ * C * (-1/3).*c1;
                                    Ax02 = eleArea * detJ * C * 1/3.*c1;
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(1,:),trtb(1,:)],...
                                                   [Ax01,Ax01,Ax01],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(2,:),trtb(2,:),trtb(2,:)],...
                                                   [Ax02,Ax02,Ax02],tsnb,trnb);
                                    Ay01 = eleArea * detJ * C * (-1/3).*c2;
                                    Ay02 = eleArea * detJ * C * 1/3.*c2;
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(1,:),trtb(1,:)],...
                                                   [Ay01,Ay01,Ay01],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(3,:),trtb(3,:),trtb(3,:)],...
                                                   [Ay02,Ay02,Ay02],tsnb,trnb);
                                case '0100'
                                    Ax01 = eleArea * detJ * C * (-1/3).*ct1;
                                    Ax02 = eleArea * detJ * C * 1/3.*ct1;
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(1,:),trtb(1,:)],...
                                                   [Ax01,Ax01,Ax01],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(2,:),trtb(2,:),trtb(2,:)],...
                                                   [Ax02,Ax02,Ax02],tsnb,trnb);
                                    Ay01 = eleArea * detJ * C * (-1/3).*ct2;
                                    Ay02 = eleArea * detJ * C * 1/3.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(1,:),trtb(1,:)],...
                                                   [Ay01,Ay01,Ay01],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],...
                                                   [trtb(3,:),trtb(3,:),trtb(3,:)],...
                                                   [Ay02,Ay02,Ay02],tsnb,trnb);
                                case '0010'
                                    A0x1 = eleArea * detJ * C * (-1/3).*c1;
                                    A0x2 = eleArea * detJ * C * 1/3.*c1;
                                    A = A + sparse([tstb(1,:),tstb(1,:),tstb(1,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0x1,A0x1,A0x1],tsnb,trnb);
                                    A = A + sparse([tstb(2,:),tstb(2,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0x2,A0x2,A0x2],tsnb,trnb);
                                    A0y1 = eleArea * detJ * C * (-1/3).*c2;
                                    A0y2 = eleArea * detJ * C * 1/3.*c2;
                                    A = A + sparse([tstb(1,:),tstb(1,:),tstb(1,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0y1,A0y1,A0y1],tsnb,trnb);
                                    A = A + sparse([tstb(3,:),tstb(3,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0y2,A0y2,A0y2],tsnb,trnb);
                                case '0001'
                                    A0x1 = eleArea * detJ * C * (-1/3).*ct1;
                                    A0x2 = eleArea * detJ * C * 1/3.*ct1;
                                    A = A + sparse([tstb(1,:),tstb(1,:),tstb(1,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0x1,A0x1,A0x1],tsnb,trnb);
                                    A = A + sparse([tstb(2,:),tstb(2,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0x2,A0x2,A0x2],tsnb,trnb);
                                    A0y1 = eleArea * detJ * C * (-1/3).*ct2;
                                    A0y2 = eleArea * detJ * C * 1/3.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(1,:),tstb(1,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0y1,A0y1,A0y1],tsnb,trnb);
                                    A = A + sparse([tstb(3,:),tstb(3,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:),trtb(3,:)],...
                                                   [A0y2,A0y2,A0y2],tsnb,trnb);
                                case '1010'
                                    Axx1 = eleArea * detJ * C * 1/2.*c1.^2;
                                    Axx2 = eleArea * detJ * C * (-1).*c1.^2;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axx1,Axx1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(2,:),Axx2,tsnb,trnb);
                                    Ayy1 = eleArea * detJ * C * 1/2.*c2.^2;
                                    Ayy2 = eleArea * detJ * C * (-1).*c2.^2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayy1,Ayy1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(3,:),Ayy2,tsnb,trnb);
                                    A = A + A';
                                    Axy1 = eleArea * detJ * C * 1.*c1.*c2;
                                    Axy2 = eleArea * detJ * C * (-1).*c1.*c2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axy1,Axy1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(2,:),trtb(1,:)],...
                                                   [Axy2,Axy2],tsnb,trnb);
                                    Ayx1 = eleArea * detJ * C * 1.*c1.*c2;
                                    Ayx2 = eleArea * detJ * C * (-1).*c1.*c2;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayx1,Ayx1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(3,:),trtb(1,:)],...
                                                   [Ayx2,Ayx2],tsnb,trnb);
                                case '1001'
                                    Axx1 = eleArea * detJ * C * 1/2.*c1.*ct1;
                                    Axx2 = eleArea * detJ * C * (-1).*c1.*ct1;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axx1,Axx1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(2,:),Axx2,tsnb,trnb);
                                    Ayy1 = eleArea * detJ * C * 1/2.*c2.*ct2;
                                    Ayy2 = eleArea * detJ * C * (-1).*c2.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayy1,Ayy1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(3,:),Ayy2,tsnb,trnb);
                                    A = A + A';
                                    Axy1 = eleArea * detJ * C * 1.*c1.*ct2;
                                    Axy2 = eleArea * detJ * C * (-1).*c1.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axy1,Axy1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(2,:),trtb(1,:)],...
                                                   [Axy2,Axy2],tsnb,trnb);
                                    Ayx1 = eleArea * detJ * C * 1.*ct1.*c2;
                                    Ayx2 = eleArea * detJ * C * (-1).*ct1.*c2;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayx1,Ayx1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(3,:),trtb(1,:)],...
                                                   [Ayx2,Ayx2],tsnb,trnb);
                                case '0110'
                                    Axx1 = eleArea * detJ * C * 1/2.*ct1.*c1;
                                    Axx2 = eleArea * detJ * C * (-1).*ct1.*c1;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axx1,Axx1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(2,:),Axx2,tsnb,trnb);
                                    Ayy1 = eleArea * detJ * C * 1/2.*ct2.*c2;
                                    Ayy2 = eleArea * detJ * C * (-1).*ct2.*c2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayy1,Ayy1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(3,:),Ayy2,tsnb,trnb);
                                    A = A + A';
                                    Axy1 = eleArea * detJ * C * 1.*ct1.*c2;
                                    Axy2 = eleArea * detJ * C * (-1).*ct1.*c2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axy1,Axy1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(2,:),trtb(1,:)],...
                                                   [Axy2,Axy2],tsnb,trnb);
                                    Ayx1 = eleArea * detJ * C * 1.*c1.*ct2;
                                    Ayx2 = eleArea * detJ * C * (-1).*c1.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayx1,Ayx1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(3,:),trtb(1,:)],...
                                                   [Ayx2,Ayx2],tsnb,trnb);
                                case '0101'
                                    Axx1 = eleArea * detJ * C * 1/2.*ct1.^2;
                                    Axx2 = eleArea * detJ * C * (-1).*ct1.^2;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axx1,Axx1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(2,:),Axx2,tsnb,trnb);
                                    Ayy1 = eleArea * detJ * C * 1/2.*ct2.^2;
                                    Ayy2 = eleArea * detJ * C * (-1).*ct2.^2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayy1,Ayy1],tsnb,trnb);
                                    A = A + sparse(tstb(1,:),trtb(3,:),Ayy2,tsnb,trnb);
                                    A = A + A';
                                    Axy1 = eleArea * detJ * C * 1.*ct1.*ct2;
                                    Axy2 = eleArea * detJ * C * (-1).*ct1.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(1,:),trtb(2,:)],...
                                                   [Axy1,Axy1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(3,:)],...
                                                   [trtb(2,:),trtb(1,:)],...
                                                   [Axy2,Axy2],tsnb,trnb);
                                    Ayx1 = eleArea * detJ * C * 1.*ct1.*ct2;
                                    Ayx2 = eleArea * detJ * C * (-1).*ct1.*ct2;
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(1,:),trtb(3,:)],...
                                                   [Ayx1,Ayx1],tsnb,trnb);
                                    A = A + sparse([tstb(1,:),tstb(2,:)],...
                                                   [trtb(3,:),trtb(1,:)],...
                                                   [Ayx2,Ayx2],tsnb,trnb);
                            end
                        case 2  % Lagrange quadratic basis
                            switch stpq
                                case '0000'
                                    
                                case '1000'

                                case '0100'

                                case '0010'

                                case '0001'

                                case '1010'

                                case '1001'

                                case '0110'

                                case '0101'

                                case '2000'

                                case '0200'

                                case '0020'

                                case '0002'

                                case '1020'

                                case '1011'

                                case '1002'

                                case '0120'

                                case '0111'

                                case '0102'

                                case '2010'

                                case '2001'

                                case '1110'

                                case '1101'

                                case '0210'

                                case '0201'

                                case '2020'

                                case '2002'

                                case '0220'

                                case '0202'

                            end
                        case 3
                            switch stpq
                                case '0000'

                                case '1000'

                                case '0100'

                                case '0010'

                                case '0001'

                                case '1010'

                                case '1001'

                                case '0110'

                                case '0101'

                            end
                    end
                case 'h'
                    switch basis.test.order
                        case 3

                        case 4

                    end
            end


        end


    case 'ipdg' % IPDG method

end
%clear dtrdx dtrdy dtsdx dtsdy ipt1 ipt2 ipt3 Aij

end
