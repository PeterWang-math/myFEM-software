function mesh = genmesh1(a,b,mesh_type,mesh_N)
%GENMESH1 此处显示有关此函数的摘要
%   此处显示详细说明
mesh.left = a; 
mesh.right = b;
mesh.type = mesh_type;
if mesh.type==3       % user-input mode
    mesh.N = mesh_N; 
    mesh.P=input('input matrix P '); 
    mesh.T=input('input matrix T '); % P & T are given by the user
elseif mesh.type==2   % adaptive FEM mode
    mesh.N = 100;
    [mesh.P,mesh.T] = genPT1(mesh);
elseif mesh.type==1   % uniform mesh mode
    mesh.N = mesh_N; 
    [mesh.P,mesh.T] = genPT1(mesh); % generate P & T matrix automatically
else                  % default mode: uniform mesh
    mesh.N = 1e3;
    [mesh.P,mesh.T] = genPT1(mesh);
end
mesh.N = size(mesh.T,2);
mesh.Nm = size(mesh.P,2); 

end

