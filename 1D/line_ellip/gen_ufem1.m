function y = gen_ufem1(mesh,basis,uh,x,femtype)
% generate FEM solution function
mesh_N = mesh.N;
mesh_P = mesh.P;
mesh_T = mesh.T;
trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;
switch femtype.type
    case 'pgm'
        for n=1:mesh_N
            if x>=mesh_P(:,mesh_T(1,n)) && x<mesh_P(:,mesh_T(2,n))
                h = mesh_P(:,mesh_T(2,n))-mesh_P(:,mesh_T(1,n));
                w0 = mesh_P(:,mesh_T(1,n))+mesh_P(:,mesh_T(2,n));
                sum = 0;
                for i=1:trnlb
                    int0 = str2func(func2str(trf{1,i}));
                    sum = sum + uh(trtb(i,n))*int0((2*x-w0)/h);
                end
                y = sum;
                break;
            elseif x == mesh_P(:,end)
                h = mesh_P(:,mesh_T(2,end))-mesh_P(:,mesh_T(1,end));
                w0 = mesh_P(:,mesh_T(1,end))+mesh_P(:,mesh_T(2,end));
                sum = 0;
                for i=1:trnlb
                    int0 = str2func(func2str(trf{1,i}));
                    sum = sum + uh(trtb(i,mesh_N))* int0((2*x-w0)/h);
                end
                y = sum;
                break;
            end
        end

    case 'ipdg'
        for n=1:mesh_N
            if n>1 && n<=mesh_N && (x==mesh_P(:,n))
                x = x-1e-15;
            end
            if x>mesh_P(:,mesh_T(1,n)) && x<mesh_P(:,mesh_T(2,n))
                h = mesh_P(:,mesh_T(2,n))-mesh_P(:,mesh_T(1,n));
                w0 = mesh_P(:,mesh_T(1,n))+mesh_P(:,mesh_T(2,n));
                sum = 0;
                for i=1:trnlb
                    int0 = str2func(func2str(trf{1,i}));
                    sum = sum + uh(trtb(i,n))*int0((2*x-w0)/h);
                end
                y = sum;
                break;
            elseif x == mesh_P(:,end)
                h = mesh_P(:,mesh_T(2,end))-mesh_P(:,mesh_T(1,end));
                w0 = mesh_P(:,mesh_T(1,end))+mesh_P(:,mesh_T(2,end));
                sum = 0;
                for i=1:trnlb
                    int0 = str2func(func2str(trf{1,i}));
                    sum = sum + uh(trtb(i,mesh_N))* int0((2*x-w0)/h);
                end
                y = sum;
                break;
            elseif x == mesh_P(:,1)
                h = mesh_P(:,mesh_T(2,1))-mesh_P(:,mesh_T(1,1));
                w0 = mesh_P(:,mesh_T(1,1))+mesh_P(:,mesh_T(2,1));
                sum = 0;
                for i=1:trnlb
                    int0 = str2func(func2str(trf{1,i}));
                    sum = sum + uh(trtb(i,1))* int0((2*x-w0)/h);
                end
                y = sum;
                break;
            end
        end
        
end


end
