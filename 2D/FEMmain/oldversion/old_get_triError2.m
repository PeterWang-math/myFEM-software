function y = get_triError2(mesh,basis,uclass,uh,type,Gauss_order,femtype)
% get error between exact solution u and FEM solution uh in triangluar mesh
% uclass : u{i}(x) = d^{i-1}u/dx^{i-1}(x) 
% uh : FEM solution, i.e. the coefficient of each global basis function
% type: 'linf' for L^infty norm error; 'l2' for L^2 norm error
%       'h1' for H^1 norm error; 'h1semi' for H^1 semi-norm error

mesh.N = size(mesh.T,2);
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
        for i=1:mesh.N
            x = linspace(mesh.P(:,i),mesh.P(:,i+1),100);
            h = mesh.P(:,mesh.T(2,i))-mesh.P(:,mesh.T(1,i));
            w0 = mesh.P(:,mesh.T(1,i))+mesh.P(:,mesh.T(2,i));
            for j=1:100
                w = 0;
                for k=1:trnlb
                    w = w + uh(trtb(k,i))*tr{1,k}((2*x(j)-w0)/h);
                end
                t = abs(feval(uclass{1,1},x(j))-w);
                y = max([y,t]);
            end
        end
    case 'l2'
        y = semi_error(mesh,basis,uclass{1,1},uh,0,0,Gauss_order);
    case 'h1semi'
        y = semi_error(mesh,basis,uclass{2},uh,1,Gauss_order);
    case 'h1'
        y1 = semi_error(mesh,basis,uclass{1},uh,0,Gauss_order);
        y2 = semi_error(mesh,basis,uclass{2},uh,1,Gauss_order);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h2semi'
        y = semi_error(mesh,basis,uclass{3},uh,2,Gauss_order);    
    case 'h2'
        y1 = get_error1(mesh,basis,uclass,uh,'h1',Gauss_order,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h2semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h3semi'
        y = semi_error(mesh,basis,uclass{4},uh,3,Gauss_order);
    case 'h3'
        y1 = get_error1(mesh,basis,uclass,uh,'h2',Gauss_order,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h3semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h4semi'
        y = semi_error(mesh,basis,uclass{5},uh,4,Gauss_order);
    case 'h4'
        y1 = get_error1(mesh,basis,uclass,uh,'h3',Gauss_order,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h4semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h5semi'
        y = semi_error(mesh,basis,uclass{6},uh,5,Gauss_order);
    case 'h5'
        y1 = get_error1(mesh,basis,uclass,uh,'h4',Gauss_order,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h5semi',Gauss_order,femtype);
        y = abs( sqrt( y1^2 + y2^2 ) );
    case 'h6semi'
        y = semi_error(mesh,basis,uclass{7},uh,6,Gauss_order);
    case 'h6'
        y1 = get_error1(mesh,basis,uclass,uh,'h5',Gauss_order,femtype);
        y2 = get_error1(mesh,basis,uclass,uh,'h6semi',Gauss_order,femtype);
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
% input d^(s+t)u/dx^sdy^t and FEM solution uh, output |u - uh|_{s+t,x/y}
y = 0;
mesh.P = mesh.P;
mesh.T = mesh.T;
mesh.N = size(mesh.T,2);
trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;
[intpt,weight] = Gauss_int_tri_ref2(Gauss_order); 
% parallel computing
for n=1:mesh.N
   J = [mesh.P(:,mesh.T(2,n))-mesh.P(:,mesh.T(1,n)), ...
        mesh.P(:,mesh.T(3,n))-mesh.P(:,mesh.T(1,n))];
   w = mesh.P(:,mesh.T(1,n));
   detJ = det(J);
   xhat = intpt*J'+w';
   f = cell(trnlb,1); 
   for k = 1:trnlb
       f{k} = @(x,y) uh(trtb(k,n)) * dphi(trf,x,y,s,t,mesh,n,k);
   end
   ipt = ( du(xhat(:,1),xhat(:,2)) - sum_fun(f,intpt(:,1),intpt(:,2)) ).^2;  
   % Rmk: in the second term above it should be sum_fun(f,intpt(:,1),intpt(:,2)),
   % not sum_fun(f,xhat(:,1),xhat(:,2)) !!!!!!
   y = y + detJ/2 * weight * ipt;  % area(ref_angle)=1/2;
end
y = abs(sqrt(y));
end


function result = sum_fun(f,x,y)
result = 0;
for i = 1:numel(f)
    result = result + f{i}(x,y);
end
end


function result = dphi(phih,x,y,s,t,mesh,n,k)
% compute d^(s+t)phi_{nk}(x,y)/dx^s dy^t
J = [mesh.P(:,mesh.T(2,n))-mesh.P(:,mesh.T(1,n)), ...
     mesh.P(:,mesh.T(3,n))-mesh.P(:,mesh.T(1,n))];
detJ = det(J);
if s==0
    if t==0
        result = phih{1,k,1}(x,y);
    else  % d^t phih_{nk}/ dy^t
        result = dphi(phih,x,y,s+1,t-1,mesh,n,k).*(-J(1,2)/detJ) + ...
                 dphi(phih,x,y,s,t-1,mesh,n,k).*(J(1,1)/detJ) ;
    end
else
    result = dphi(phih,x,y,s,t,mesh,n,k);
end
    
end


