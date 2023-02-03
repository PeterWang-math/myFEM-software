function F = assembQuadload(mesh,basis,f,p,q,Gauss_order,femtype)
% assemble load vector F for 2d case
% assemble int_{En} f x d^{p+q}psai/dx^{p}dy^{q}

tsnb = basis.test.Nb; 
tsnlb = basis.test.Nlb;
tsf = basis.test.reffun;
tstb = basis.test.Tb;

F = sparse(tsnb,1);
if isreal(f) && (f==0) % f = 0;
    return
end

switch size(mesh.T,1)
    case 3
        [intpt,weight] = Gauss_int_tri_ref2(Gauss_order); 
        eleArea = 1/2;
    case 4
        [intpt,weight] = Gauss_int_quad_ref2(Gauss_order); 
        eleArea = 1;
end


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

        if ~isnumeric(f)
        % f is a function, use Gauss quadrature to compute integral
            xhat = J*intpt'+reshape(pn1,[],1);
            xhat = xhat';
            pxy = reshape(xhat,[],mesh.N);
            px = reshape(pxy(1:length(weight),:),[],1);
            py = reshape(pxy(length(weight)+1:end,:),[],1);
            clear pxy xhat
            for j = 1:tsnlb
                ipt = f(px,py).*repmat(tsf{j,p+1,p+q+1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                dj = eleArea * detJ.* (weight * reshape(ipt,[],mesh.N));
                F = F + sparse(tstb(j,:),1,dj,tsnb,1);
            end
        elseif eleArea == 1/2  
        % f is a constant, use barycentric coordinates to compute integral
        % only works for triangluar mesh
            switch basis.test.type
                case 'l'
                    switch basis.test.order
                        case 1  % Lagrange linear basis
                            d1 = eleArea * detJ * f * 1/3;
                            d2 = eleArea * detJ * f * 1/3;
                            d3 = eleArea * detJ * f * 1/3;
                            F = F + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],1,[d1,d2,d3],tsnb,1);
                        case 2  % Lagrange quadratic basis
                            d1 = eleArea * detJ * f * 0;
                            d2 = eleArea * detJ * f * 1/3;
                            d3 = eleArea * detJ * f * 1/3;
                            d4 = eleArea * detJ * f * 0;
                            d5 = eleArea * detJ * f * 1/3;
                            d6 = eleArea * detJ * f * 0;
                            F = F + sparse([tstb(1,:),tstb(2,:),tstb(3,:),tstb(4,:),tstb(5,:),tstb(6,:)],...
                                1,[d1,d2,d3,d4,d5,d6],tsnb,1);
                        case 3

                    end
                case 'h'    
                    switch basis.test.order
                        case 3

                        case 4

                    end
            end
        elseif eleArea == 1
        % f is a constant, use barycentric coordinates to compute integral
        % only works for quadrilateral mesh
            switch basis.test.type
                case 'l'
                    switch basis.test.order
                        case 1  % Lagrange linear basis
                            
                        case 2  % Lagrange quadratic basis
                            
                        case 3

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
F = full(F);
% clear px py ipt dj;

end
