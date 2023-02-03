function [P,T] = genPT1(mesh)
% generate infomation matrix P,T of uniform 1d mesh 
% input info: a struct 'mesh' consisting mesh.left, mesh.right and mesh.N
% mesh.left: the left node of the interval, i.e., a 
% mesh.right: the right node of the interval, i.e. b
% mesh.N: number of mesh elements

a = mesh.left; b = mesh.right; N = mesh.N;
P = a:(b-a)/N:b;
T = zeros(2,N);
T(1,:)=1:1:N;
T(2,:)=2:1:N+1;

end

