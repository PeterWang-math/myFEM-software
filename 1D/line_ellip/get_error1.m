function y = get_error1(mesh,basis,uclass,uh,type,Gauss_nodes,femtype)
% plot error between exact solution u and FEM solution uh
% uclass : u{i}(x) = d^{i-1}u/dx^{i-1}(x) 
% uh : FEM solution, i.e. the coefficient of each global basis function
% type: 'linf' for L^infty norm error; 'l2' for L^2 norm error
%       'h1' for H^1 norm error; 'h1semi' for H^1 semi-norm error
mesh_N = mesh.N;
mesh_P = mesh.P;
mesh_T = mesh.T;
trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;
tr = cell(1,trnlb);
for k=1:trnlb
    tr{1,k} = str2func(func2str(trf{1,k}));
end

switch type
    case 'linf' % select 100 points on each element to estimate error
        y = 0; 
        for i=1:mesh_N
            x = linspace(mesh_P(:,i),mesh_P(:,i+1),100);
            h = mesh_P(:,mesh_T(2,i))-mesh_P(:,mesh_T(1,i));
            w0 = mesh_P(:,mesh_T(1,i))+mesh_P(:,mesh_T(2,i));
            for j=1:100
                w = 0;
                for k=1:trnlb
                    w = w + uh(trtb(k,i))*tr{1,k}((2*x(j)-w0)/h);
                end
                t = abs(feval(uclass{1},x(j))-w);
                y = max([y,t]);
            end
        end
    case 'l2'
        y = semi_error(mesh,basis,uclass{1},uh,0,Gauss_nodes);
    case 'h1semi'
        y = semi_error(mesh,basis,uclass{2},uh,1,Gauss_nodes);
    case 'h1'
        y1 = semi_error(mesh,basis,uclass{1},uh,0,Gauss_nodes);
        y2 = semi_error(mesh,basis,uclass{2},uh,1,Gauss_nodes);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h2semi'
        y = semi_error(mesh,basis,uclass{3},uh,2,Gauss_nodes);    
    case 'h2'
        y1 = get_error1(mesh,basis,uclass,uh,'h1',Gauss_nodes,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h2semi',Gauss_nodes,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h3semi'
        y = semi_error(mesh,basis,uclass{4},uh,3,Gauss_nodes);
    case 'h3'
        y1 = get_error1(mesh,basis,uclass,uh,'h2',Gauss_nodes,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h3semi',Gauss_nodes,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h4semi'
        y = semi_error(mesh,basis,uclass{5},uh,4,Gauss_nodes);
    case 'h4'
        y1 = get_error1(mesh,basis,uclass,uh,'h3',Gauss_nodes,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h4semi',Gauss_nodes,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h5semi'
        y = semi_error(mesh,basis,uclass{6},uh,5,Gauss_nodes);
    case 'h5'
        y1 = get_error1(mesh,basis,uclass,uh,'h4',Gauss_nodes,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h5semi',Gauss_nodes,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h6semi'
        y = semi_error(mesh,basis,uclass{7},uh,6,Gauss_nodes);
    case 'h6'
        y1 = get_error1(mesh,basis,uclass,uh,'h5',Gauss_nodes,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h6semi',Gauss_nodes,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
end

% if strcmp(femtype.type,'ipdg') == 1
%     switch type
%         case 'ipdg_h1'
%             y = semi_error(mesh,basis,uclass{1},uh,0,Gauss_nodes);
%             
%     end
% end


end



function y = semi_error(mesh,basis,du,uh,s,Gauss_nodes)
% input d^su/dx^s and FEM solution uh, output |u - uh|_{H^s}
y = 0;
mesh_N = mesh.N;
mesh_P = mesh.P;
mesh_T = mesh.T;
trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;
[intpt,weight] = Gauss_int_ref1(Gauss_nodes); 
% parallel computing
for n=1:mesh_N
   h = mesh_P(:,mesh_T(2,n))-mesh_P(:,mesh_T(1,n));
   w0 = mesh_P(:,mesh_T(1,n))+mesh_P(:,mesh_T(2,n));
   w = w0/2;
   f = cell(trnlb,1); 
   for k = 1:trnlb
       int0 = str2func(func2str(trf{s+1,k}));
       f{k} = @(x) uh(trtb(k,n))*int0((2*x-w0)/h);
   end 
%    vec = zeros(1,Gauss_nodes);
%    for k=1:Gauss_nodes
%        pt = int_point(k);
%        vec(k) = ( du(h/2*pt+w) - (2/h)^s * sum_fun(f,h/2*pt+w) )^2;
%    end
   ipt = ( du(h/2*intpt+w) - (2/h)^s * sum_fun(f,h/2*intpt+w) ).^2;
   y = y + h/2 * weight * ipt';
end
y = abs(sqrt(y));
end

function y = sum_fun(f,x)
y = 0;
for i = 1:numel(f)
    y = y + f{i}(x);
end
end

