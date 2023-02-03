function F = assembTriload(mesh,basis,f,p,q,Gauss_order,femtype)
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

pq = strcat([num2str(p),num2str(q)]);
[intpt,weight] = Gauss_int_tri_ref2(Gauss_order); 
eleArea = 1/2;

switch femtype.type
    case 'pgm' % PG method
        pn1 = mesh.P(:,mesh.T(1,:));
        v21 = mesh.P(:,mesh.T(2,:))-pn1;
        v31 = mesh.P(:,mesh.T(3,:))-pn1;
        detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);
        J = [reshape(v21,[],1),reshape(v31,[],1)];

        if ~isnumeric(f)  % f is a function, use Gauss quadrature to compute integral
            xhat = J*intpt'+reshape(pn1,[],1);
            xhat = xhat';
            pxy = reshape(xhat,[],mesh.N);
            px = reshape(pxy(1:length(weight),:),[],1);
            py = reshape(pxy(length(weight)+1:end,:),[],1);
            clear pxy xhat
            for j = 1:tsnlb
                switch pq
                    case '00'
                        ipt = f(px,py).*repmat(tsf{j,1,1}(intpt(:,1),intpt(:,2)),mesh.N,1);
                        dj = eleArea * detJ.* (weight * reshape(ipt,[],mesh.N));
                        F = F + sparse(tstb(j,:),1,dj,tsnb,1);
                    case '10'

                    case '01'

                    case '20'

                    case '11'

                    case '02'

                end
            end
        else  % f is a constant, use barycentric coordinates to compute integral
            switch basis.test.type
                case 'l'
                    switch basis.test.order
                        case 1  % Lagrange linear basis
                            switch pq
                                case '00'
                                    d1 = eleArea * detJ * f * 1/3;
                                    d2 = eleArea * detJ * f * 1/3;
                                    d3 = eleArea * detJ * f * 1/3;
                                    F = F + sparse([tstb(1,:),tstb(2,:),tstb(3,:)],1,...
                                        [d1,d2,d3],tsnb,1);
                                case '10'

                                case '01'

                            end
                        case 2  % Lagrange quadratic basis
                            switch pq
                                case '00'
                                    d2 = eleArea * detJ * f * 1/3;
                                    d3 = eleArea * detJ * f * 1/3;
                                    d5 = eleArea * detJ * f * 1/3;
                                    F = F + sparse([tstb(2,:),tstb(3,:),tstb(5,:)],...
                                        1,[d2,d3,d5],tsnb,1);
                                case '10'

                                case '01'

                                case '20'

                                case '11'

                                case '02'

                            end
                        case 3  % Lagrange cubic basis

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
