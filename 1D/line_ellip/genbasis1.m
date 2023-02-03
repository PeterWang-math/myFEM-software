function newbasis = genbasis1(mesh,basis,femtype)
% generate basis function space infomation and FEM nodes infomation in 1d,
% trial and test space are considered seperately for PGM
%
% The following input infomation is needed:
% mesh: a struct consisting of mesh.left, mesh.right, mesh.P, mesh.T,  
%       mesh.N, mesh.Nm;
% basis: a struct consisting of basis.trial/test.type, basis.trial/test.order
%
% The output information includes:
% basis.trial/test.reffun: a 7 x Nlb cell consisting of all the local basis 
%                        functions on reference element, the ith row stores the
%                        i-1 derivative of the basis function, while
%					     the column is corresponding to the nodal basis                      
% basis.trial/test.Pb/Tb/Nb/Nlb: FEM nodes information

% linear nodal basis 
%l1 = @(x) 1/2*(-x+1);
%l2 = @(x) 1/2*(x+1);

% generate trial basis function as well as Pb and Tb matrix, 
% the basis function is defined on the reference element (-1,1)
if strcmp(femtype.type,'pgm')==1 || strcmp(femtype.type,'ipdg')==1 % PG/IPDG method
   
if basis.trial.type ~= 'l' && basis.trial.type ~= 'h'
    warning('The basis.trial.type must be either l or h, already reset to l');
    basis.trial.type = 'l';
end
if basis.trial.type == 'l'
    ord = basis.trial.order;
    if floor(ord)~=ord || ord<1 || ord>6
        warning(['The basis.trial.order must be a positive integer less than 7,',...
                'already reset to 1']);
        basis.trial.order = 1;
    end
    switch basis.trial.order
          case 1 % Lagrange linear basis function
              switch femtype.type
                  case 'pgm'
                      basis.trial.Pb = mesh.P;
                      basis.trial.Tb = mesh.T;
                  case 'ipdg'
                      Pb = mesh.P;
                      Pb2 = Pb(:,2:end-1);
                      basis.trial.Pb = [Pb(:,1),reshape([Pb2;Pb2],1,2*mesh.N-2),Pb(:,end)];
                      basis.trial.Tb = reshape(1:2*mesh.N,2,mesh.N);
              end
                 fclass = cell(7,2);
%                  fclass{1,1} = @(x) l1(x);
%                  fclass{1,2} = @(x) l2(x);
                 fclass{1,1} = @(x) 1/2*(-x+1);
                 fclass{1,2} = @(x) 1/2*(x+1);
                 fclass{2,1} = @(x) -1/2;
                 fclass{2,2} = @(x) 1/2;
                 fclass{3,1} = @(x) 0;  fclass{3,2} = @(x) 0;
                 fclass{4,1} = @(x) 0;  fclass{4,2} = @(x) 0;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;
                 basis.trial.reffun = fclass; 
          case 2 % Lagrange quadratic basis function
              switch femtype.type
                  case 'pgm'
                      basis.trial.Pb = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                      Tb = zeros(3,mesh.N);
                      Tb(1,:) = 1:2:2*mesh.N-1; Tb(2,:) = 2:2:2*mesh.N;
                      Tb(3,:) = 3:2:2*mesh.N+1;
                      basis.trial.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                      Pb = zeros(3,mesh.N);
                      Pb(1,:) = 1:2:2*mesh.N-1; Pb(2,:) = 2:2:2*mesh.N;
                      Pb(3,:) = 3:2:2*mesh.N+1;
                      Pb = reshape(Pb,1,3*mesh.N);
                      basis.trial.Pb = Pb0(Pb);
                      basis.trial.Tb = reshape(1:3*mesh.N,3,mesh.N);
               end
                 fclass = cell(7,3);
%                  fclass{1,1} = @(x) l1(x)*(l1(x)-l2(x));
%                  fclass{1,2} = @(x) 4*l1(x)*l2(x);
%                  fclass{1,3} = @(x) l2(x)*(l2(x)-l1(x));
                 fclass{1,1} = @(x) 1/2*(-1 + x).*x;
                 fclass{1,2} = @(x) 1-x.^2;
                 fclass{1,3} = @(x) 1/2*x.*(1 + x);
                 fclass{2,1} = @(x) -1/2 + x;
                 fclass{2,2} = @(x) -2*x;
                 fclass{2,3} = @(x) 1/2 + x;
                 fclass{3,1} = @(x) 1;
                 fclass{3,2} = @(x) -2;
                 fclass{3,3} = @(x) 1;
                 fclass{4,1} = @(x) 0;  fclass{4,2} = @(x) 0;  fclass{4,3} = @(x) 0;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;  fclass{5,3} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  fclass{6,3} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  fclass{7,3} = @(x) 0;
                 basis.trial.reffun = fclass;
          case 3 % Lagrange cubic basis function
              switch femtype.type
                  case 'pgm'
                      basis.trial.Pb = mesh.left:(mesh.right-mesh.left)/(3*mesh.N):mesh.right;
                      Tb = zeros(4,mesh.N);
                      Tb(1,:) = 1:3:3*mesh.N-2; Tb(2,:) = 2:3:3*mesh.N-1; 
                      Tb(3,:) = 3:3:3*mesh.N; Tb(4,:) = 4:3:3*mesh.N+1;
                      basis.trial.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(3*mesh.N):mesh.right;
                      Pb = zeros(4,mesh.N);
                      Pb(1,:) = 1:3:3*mesh.N-2; Pb(2,:) = 2:3:3*mesh.N-1; 
                      Pb(3,:) = 3:3:3*mesh.N; Pb(4,:) = 4:3:3*mesh.N+1;
                      Pb = reshape(Pb,1,4*mesh.N);
                      basis.trial.Pb = Pb0(Pb);
                      basis.trial.Tb = reshape(1:4*mesh.N,4,mesh.N);
              end
                 fclass = cell(7,4);
%                  fclass{1,1} = @(x) 1/2*l1(x)*(l1(x)-2*l2(x))*(2*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 9/2*l1(x)*l2(x)*(2*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -9/2*l1(x)*l2(x)*(l1(x)-2*l2(x));
%                  fclass{1,4} = @(x) 1/2*l2(x)*(l1(x)-2*l2(x))*(2*l1(x)-l2(x));
                 fclass{1,1} = @(x) -1/16*(-1 + x).*(-1 + 3*x).*(1 + 3*x);
                 fclass{1,2} = @(x) 9/16*(-1 + x).*(1 + x).*(-1 + 3*x);
                 fclass{1,3} = @(x) -9/16*(-1 + x).*(1 + x).*(1 + 3*x);
                 fclass{1,4} = @(x) 1/16*(1 + x).*(-1 + 3*x).*(1 + 3*x);
                 fclass{2,1} = @(x) 1/16*(1 + 18*x - 27*x.^2);
                 fclass{2,2} = @(x) 9/16*(-3 - 2*x + 9*x.^2);
                 fclass{2,3} = @(x) -9/16*(-3 + 2*x + 9*x.^2);
                 fclass{2,4} = @(x) 1/16*(-1 + 18*x + 27*x.^2);
                 fclass{3,1} = @(x) -9/8*(-1 + 3*x);
                 fclass{3,2} = @(x) 9/8*(-1 + 9*x);
                 fclass{3,3} = @(x) -9/8*(1 + 9*x);
                 fclass{3,4} = @(x) 9/8*(1 + 3*x);
                 fclass{4,1} = @(x) -27/8;  
                 fclass{4,2} = @(x) 81/8;  
                 fclass{4,3} = @(x) -81/8;
                 fclass{4,4} = @(x) 27/8;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;  
                 fclass{5,3} = @(x) 0;  fclass{5,4} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 basis.trial.reffun = fclass;
          case 4 % Lagrange fourth-order basis function
              switch femtype.type
                  case 'pgm'
                      basis.trial.Pb = mesh.left:(mesh.right-mesh.left)/(4*mesh.N):mesh.right;
                      Tb = zeros(5,mesh.N);
                      Tb(1,:) = 1:4:4*mesh.N-3; Tb(2,:) = 2:4:4*mesh.N-2; 
                      Tb(3,:) = 3:4:4*mesh.N-1; Tb(4,:) = 4:4:4*mesh.N;
                      Tb(5,:) = 5:4:4*mesh.N+1;
                      basis.trial.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(4*mesh.N):mesh.right;
                      Pb = zeros(5,mesh.N);
                      Pb(1,:) = 1:4:4*mesh.N-3; Pb(2,:) = 2:4:4*mesh.N-2; 
                      Pb(3,:) = 3:4:4*mesh.N-1; Pb(4,:) = 4:4:4*mesh.N;
                      Pb(5,:) = 5:4:4*mesh.N+1;
                      Pb = reshape(Pb,1,5*mesh.N);
                      basis.trial.Pb = Pb0(Pb);
                      basis.trial.Tb = reshape(1:5*mesh.N,5,mesh.N);
              end
                 fclass = cell(7,5);
%                  fclass{1,1} = @(x) 1/3*l1(x)*(l1(x)-l2(x))*(l1(x)-3*l2(x))*(3*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 16/3*l1(x)*l2(x)*(l1(x)-l2(x))*(3*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -4*l1(x)*l2(x)*(l1(x)-3*l2(x))*(3*l1(x)-l2(x));
%                  fclass{1,4} = @(x) 16/3*l1(x)*l2(x)*(l1(x)-l2(x))*(l1(x)-3*l2(x));
%                  fclass{1,5} = @(x) -1/3*l2(x)*(l1(x)-l2(x))*(l1(x)-3*l2(x))*(3*l1(x)-l2(x));
                 fclass{1,1} = @(x) 1/6*(-1 + x).*x.*(-1 + 2*x).*(1 + 2*x);
                 fclass{1,2} = @(x) -4/3*(-1 + x).*x.*(1 + x).*(-1 + 2*x);
                 fclass{1,3} = @(x) 1 - 5*x.^2 + 4*x.^4;
                 fclass{1,4} = @(x) -4/3*(-1 + x).*x.*(1 + x).*(1 + 2*x);
                 fclass{1,5} = @(x) 1/6*x.*(1 + x).*(-1 + 2*x).*(1 + 2*x);
                 fclass{2,1} = @(x) 1/6*(1 - 2*x - 12*x.^2 + 16*x.^3);
                 fclass{2,2} = @(x) -4/3*(1 - 4*x - 3*x.^2 + 8*x.^3);
                 fclass{2,3} = @(x) 2*x.*(-5 + 8*x.^2);
                 fclass{2,4} = @(x) -4/3*(-1 - 4*x + 3*x.^2 + 8*x.^3);
                 fclass{2,5} = @(x) 1/6*(-1 - 2*x + 12*x.^2 + 16*x.^3);
                 fclass{3,1} = @(x) -1/3 - 4*x + 8*x.^2;
                 fclass{3,2} = @(x) 16/3 + 8*x - 32*x.^2;
                 fclass{3,3} = @(x) -10 + 48*x.^2;
                 fclass{3,4} = @(x) 16/3 - 8*x - 32*x.^2;
                 fclass{3,5} = @(x) -1/3 + 4*x + 8*x.^2;
                 fclass{4,1} = @(x) -4 + 16*x;  
                 fclass{4,2} = @(x) 8 - 64*x;  
                 fclass{4,3} = @(x) 96*x;
                 fclass{4,4} = @(x) -8*(1 + 8*x);
                 fclass{4,5} = @(x) 4 + 16*x;
                 fclass{5,1} = @(x) 16;  
                 fclass{5,2} = @(x) -64;  
                 fclass{5,3} = @(x) 96;  
                 fclass{5,4} = @(x) -64; 
                 fclass{5,5} = @(x) 16;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{6,5} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;
                 basis.trial.reffun = fclass;
          case 5 % Lagrange fifth-order basis function
              switch femtype.type
                  case 'pgm'
                      basis.trial.Pb = mesh.left:(mesh.right-mesh.left)/(5*mesh.N):mesh.right;
                      Tb = zeros(6,mesh.N);
                      Tb(1,:) = 1:5:5*mesh.N-4; Tb(2,:) = 2:5:5*mesh.N-3; 
                      Tb(3,:) = 3:5:5*mesh.N-2; Tb(4,:) = 4:5:5*mesh.N-1;
                      Tb(5,:) = 5:5:5*mesh.N; Tb(6,:) = 6:5:5*mesh.N+1;
                      basis.trial.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(5*mesh.N):mesh.right;
                      Pb = zeros(6,mesh.N);
                      Pb(1,:) = 1:5:5*mesh.N-4; Pb(2,:) = 2:5:5*mesh.N-3; 
                      Pb(3,:) = 3:5:5*mesh.N-2; Pb(4,:) = 4:5:5*mesh.N-1;
                      Pb(5,:) = 5:5:5*mesh.N; Pb(6,:) = 6:5:5*mesh.N+1;
                      Pb = reshape(Pb,1,6*mesh.N);
                      basis.trial.Pb = Pb0(Pb);
                      basis.trial.Tb = reshape(1:6*mesh.N,6,mesh.N);
              end
                 fclass = cell(7,6);
%                  fclass{1,1} = @(x) 1/24*l1(x)*(l1(x)-4*l2(x))*(2*l1(x)-3*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 25/24*l1(x)*l2(x)*(2*l1(x)-3*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -25/12*l1(x)*l2(x)*(l1(x)-4*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,4} = @(x) 25/12*l1(x)*l2(x)*(l1(x)-4*l2(x))...
%                                     *(2*l1(x)-3*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,5} = @(x) -25/24*l1(x)*l2(x)*(l1(x)-4*l2(x))...
%                                     *(2*l1(x)-3*l2(x))*(3*l1(x)-2*l2(x));
%                  fclass{1,6} = @(x) 1/24*l2(x)*(l1(x)-4*l2(x))*(2*l1(x)-3*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
                 fclass{1,1} = @(x) -1/768*(-1 + x).*(-3 + 5*x).*(-1 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{1,2} = @(x) 25/768*(-1 + x).*(1 + x).*(-3 + 5*x).*(-1 + 5*x).*(1 + 5*x);
                 fclass{1,3} = @(x) -25/384*(-1 + x).*(1 + x).*(-3 + 5*x).*(-1 + 5*x).*(3 + 5*x);
                 fclass{1,4} = @(x) 25/384*(-1 + x).*(1 + x).*(-3 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{1,5} = @(x) -25/768*(-1 + x).*(1 + x).*(-1 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{1,6} = @(x) 1/768*(1 + x).*(-3 + 5*x).*(-1 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{2,1} = @(x) 1/768*(-9 - 500*x + 750*x.^2 + 2500*x.^3 - 3125*x.^4);
                 fclass{2,2} = @(x) 25/768*(5 + 156*x - 390*x.^2 - 300*x.^3 + 625*x.^4);
                 fclass{2,3} = @(x) -25/384*(45 + 68*x - 510*x.^2 - 100*x.^3 + 625*x.^4);
                 fclass{2,4} = @(x) 25/384*(45 - 68*x - 510*x.^2 + 100*x.^3 + 625*x.^4);
                 fclass{2,5} = @(x) -25/768*(5 - 156*x - 390*x.^2 + 300*x.^3 + 625*x.^4);
                 fclass{2,6} = @(x) 1/768*(9 - 500*x - 750*x.^2 + 2500*x.^3 + 3125*x.^4);
                 fclass{3,1} = @(x) -125/192*(1 - 3*x - 15*x.^2 + 25*x.^3);
                 fclass{3,2} = @(x) 25/192*(39 - 195*x - 225*x.^2 + 625*x.^3);
                 fclass{3,3} = @(x) -25/96*(17 - 255*x - 75*x.^2 + 625*x.^3);
                 fclass{3,4} = @(x) 25/96*(-17 - 255*x + 75*x.^2 + 625*x.^3);
                 fclass{3,5} = @(x) -25/192*(-39 - 195*x + 225*x.^2 + 625*x.^3);
                 fclass{3,6} = @(x) 125/192*(-1 - 3*x + 15*x.^2 + 25*x.^3);
                 fclass{4,1} = @(x) -125/64*(-1 - 10*x + 25*x.^2);  
                 fclass{4,2} = @(x) 125/64*(-13 - 30*x + 125*x.^2);  
                 fclass{4,3} = @(x) -125/32*(-17 - 10*x + 125*x.^2);
                 fclass{4,4} = @(x) 125/32*(-17 + 10*x + 125*x.^2);
                 fclass{4,5} = @(x) -125/64*(-13 + 30*x + 125*x.^2);
                 fclass{4,6} = @(x) 125/64*(-1 + 10*x + 25*x.^2);
                 fclass{5,1} = @(x) -625/32*(-1 + 5*x);  
                 fclass{5,2} = @(x) 625/32*(-3 + 25*x);  
                 fclass{5,3} = @(x) -625/16*(-1 + 25*x);  
                 fclass{5,4} = @(x) 625/16*(1 + 25*x); 
                 fclass{5,5} = @(x) -625/32*(3 + 25*x);
                 fclass{5,6} = @(x) 625/32*(1 + 5*x);
                 fclass{6,1} = @(x) -3125/32;  
                 fclass{6,2} = @(x) 15625/32;  
                 fclass{6,3} = @(x) -15625/16;  
                 fclass{6,4} = @(x) 15625/16;
                 fclass{6,5} = @(x) -15625/32;  
                 fclass{6,6} = @(x) 3125/32;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;  fclass{7,6} = @(x) 0;
                 basis.trial.reffun = fclass;
          case 6 % Lagrange sixth-order basis function
               switch femtype.type
                  case 'pgm'
                      basis.trial.Pb = mesh.left:(mesh.right-mesh.left)/(6*mesh.N):mesh.right;
                      Tb = zeros(7,mesh.N);
                      Tb(1,:) = 1:6:6*mesh.N-5; Tb(2,:) = 2:6:6*mesh.N-4; 
                      Tb(3,:) = 3:6:6*mesh.N-3; Tb(4,:) = 4:6:6*mesh.N-2;
                      Tb(5,:) = 5:6:6*mesh.N-1; Tb(6,:) = 6:6:6*mesh.N;
                      Tb(7,:) = 7:6:6*mesh.N+1;
                      basis.trial.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(6*mesh.N):mesh.right;
                      Pb = zeros(7,mesh.N);
                      Pb(1,:) = 1:6:6*mesh.N-5; Pb(2,:) = 2:6:6*mesh.N-4; 
                      Pb(3,:) = 3:6:6*mesh.N-3; Pb(4,:) = 4:6:6*mesh.N-2;
                      Pb(5,:) = 5:6:6*mesh.N-1; Pb(6,:) = 6:6:6*mesh.N;
                      Pb(7,:) = 7:6:6*mesh.N+1;
                      Pb = reshape(Pb,1,7*mesh.N);
                      basis.trial.Pb = Pb0(Pb);
                      basis.trial.Tb = reshape(1:7*mesh.N,7,mesh.N);
               end
                 fclass = cell(7,7);
%                  fclass{1,1} = @(x) 1/10*l1(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 18/5*l1(x)*l2(x)*(l1(x)-2*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -9/2*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,4} = @(x) 4*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,5} = @(x) -9/2*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))...
%                                     *(l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,6} = @(x) 18/5*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))...
%                                     *(l1(x)-l2(x))*(2*l1(x)-l2(x));
%                  fclass{1,7} = @(x) -1/10*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
                 fclass{1,1} = @(x) 1/80*(-1 + x).*x.*(-2 + 3*x).*(-1 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{1,2} = @(x) -9/40*(-1 + x).*x.*(1 + x).*(-2 + 3*x).*(-1 + 3*x).*(1 + 3*x);
                 fclass{1,3} = @(x) 9/16*(-1 + x).*x.*(1 + x).*(-2 + 3*x).*(-1 + 3*x).*(2 + 3*x);
                 fclass{1,4} = @(x) 1/4*(4 - 49*x.^2 + 126*x.^4 - 81*x.^6);
                 fclass{1,5} = @(x) 9/16*(-1 + x).*x.*(1 + x).*(-2 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{1,6} = @(x) -9/40*(-1 + x).*x.*(1 + x).*(-1 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{1,7} = @(x) 1/80*x.*(1 + x).*(-2 + 3*x).*(-1 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{2,1} = @(x) 1/80*(-4 + 8*x + 135*x.^2 - 180*x.^3 - 405*x.^4 + 486*x.^5);
                 fclass{2,2} = @(x) -9/20*(-1 + 3*x + 30*x.^2 - 60*x.^3 - 45*x.^4 + 81*x.^5);
                 fclass{2,3} = @(x) 9/16*(-4 + 24*x + 39*x.^2 - 156*x.^3 - 45*x.^4 + 162*x.^5);
                 fclass{2,4} = @(x) 1/2*(-49*x + 252*x.^3 - 243*x.^5);
                 fclass{2,5} = @(x) 9/16*(4 + 24*x - 39*x.^2 - 156*x.^3 + 45*x.^4 + 162*x.^5);
                 fclass{2,6} = @(x) -9/20*(1 + 3*x - 30*x.^2 - 60*x.^3 + 45*x.^4 + 81*x.^5);
                 fclass{2,7} = @(x) 1/80*(4 + 8*x - 135*x.^2 - 180*x.^3 + 405*x.^4 + 486*x.^5);
                 fclass{3,1} = @(x) 1/40*(4 + 135*x - 270*x.^2 - 810*x.^3 + 1215*x.^4);
                 fclass{3,2} = @(x) -27/20*(1 + 20*x - 60*x.^2 - 60*x.^3 + 135*x.^4);
                 fclass{3,3} = @(x) 27/8*(4 + 13*x - 78*x.^2 - 30*x.^3 + 135*x.^4);
                 fclass{3,4} = @(x) 1/2*(-49 + 756*x.^2 - 1215*x.^4);
                 fclass{3,5} = @(x) 27/8*(4 - 13*x - 78*x.^2 + 30*x.^3 + 135*x.^4);
                 fclass{3,6} = @(x) -27/20*(1 - 20*x - 60*x.^2 + 60*x.^3 + 135*x.^4);
                 fclass{3,7} = @(x) 1/40*(4 - 135*x - 270*x.^2 + 810*x.^3 + 1215*x.^4);
                 fclass{4,1} = @(x) 27/8*(1 - 4*x - 18*x.^2 + 36*x.^3);  
                 fclass{4,2} = @(x) -27*(1 - 6*x - 9*x.^2 + 27*x.^3);  
                 fclass{4,3} = @(x) 27/8*(13 - 156*x - 90*x.^2 + 540*x.^3);
                 fclass{4,4} = @(x) 54*x.*(14 - 45*x.^2);
                 fclass{4,5} = @(x) 27/8*(-13 - 156*x + 90*x.^2 + 540*x.^3);
                 fclass{4,6} = @(x) -27*(-1 - 6*x + 9*x.^2 + 27.*x.^3);
                 fclass{4,7} = @(x) 27/8*(-1 - 4*x + 18*x.^2 + 36*x.^3);
                 fclass{5,1} = @(x) 27/2*(-1 - 9*x + 27*x.^2);  
                 fclass{5,2} = @(x) -81*(-2 - 6*x + 27*x.^2);  
                 fclass{5,3} = @(x) 81/2*(-13 - 15*x + 135*x.^2);  
                 fclass{5,4} = @(x) 756 - 7290*x.^2; 
                 fclass{5,5} = @(x) 81/2*(-13 + 15*x + 135*x.^2);
                 fclass{5,6} = @(x) -81*(-2 + 6*x + 27*x.^2);
                 fclass{5,7} = @(x) 27/2*(-1 + 9*x + 27*x.^2);
                 fclass{6,1} = @(x) 243/2*(-1 + 6*x);  
                 fclass{6,2} = @(x) -486*(-1 + 9*x);  
                 fclass{6,3} = @(x) 1215/2*(-1 + 18*x);  
                 fclass{6,4} = @(x) -14580*x;
                 fclass{6,5} = @(x) 1215/2*(1 + 18*x);  
                 fclass{6,6} = @(x) -486*(1 + 9*x);
                 fclass{6,7} = @(x) 243/2*(1 + 6*x);
                 fclass{7,1} = @(x) 729;  
                 fclass{7,2} = @(x) -4374;  
                 fclass{7,3} = @(x) 10935;  
                 fclass{7,4} = @(x) -14580;
                 fclass{7,5} = @(x) 10935;  
                 fclass{7,6} = @(x) -4374;
                 fclass{7,7} = @(x) 729;
                 basis.trial.reffun = fclass;      
    end
    basis.trial.Nb = size(basis.trial.Pb,2);
    basis.trial.Nlb = size(basis.trial.Tb,1);
elseif basis.trial.type == 'h'
    ord = basis.trial.order;
    if ord<3 || ord>6
        warning(['The basis.trial.order must be an interger more than 2 less than 7 ,',...
                'already reset to 3']);
        basis.trial.order = 3;
    end
        switch basis.trial.order  
            case 3  % Hermite cubic element: a1, da1, a2, da2
               switch femtype.type
                  case 'pgm' % local basis index for 1st element: 1,2,N+2,N+3
                      Pb = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                      basis.trial.Pb = reshape([Pb; Pb],1,2*mesh.N+2);
                      Tb = zeros(4,mesh.N);
                      Tb(1,:) = 1:2:2*mesh.N-1; Tb(2,:) = Tb(1,:)+2; 
                      Tb(3,:) = Tb(1,:)+1; Tb(4,:) = Tb(1,:)+3; 
                      basis.trial.Tb = Tb;
                  case 'ipdg' % local basis index for 1st element: 1,2,3,4
                      Pb0 = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                      Tb = reshape(1:4*mesh.N,4,mesh.N);
                      t = Tb(2,:); Tb(2,:) = Tb(3,:); Tb(3,:) = t;
                      basis.trial.Tb = Tb;
                      Pb = zeros(4,mesh.N);
                      Pb(1,:) = 1:mesh.N;  Pb(2,:) = Pb(1,:);
                      Pb(3,:) = 2:mesh.N+1; Pb(4,:) = Pb(3,:);
                      basis.trial.Pb = Pb0(reshape(Pb,1,4*mesh.N));
               end
                 fclass = cell(7,4);
%                  fclass{1,1} = @(x) l1(x)^2*(l1(x)+3*l2(x));
%                  fclass{1,2} = @(x) l2(x)^2*(3*l1(x)+l2(x));
%                  fclass{1,3} = @(x) 2*l1(x)^2*l2(x);
%                  fclass{1,4} = @(x) -2*l1(x)*l2(x)^2;
                 fclass{1,1} = @(x) 1/4*(-1 + x).^2.*(2 + x);
                 fclass{1,2} = @(x) -1/4*(-2 + x).*(1 + x).^2;
                 fclass{1,3} = @(x) 1/4*(-1 + x).^2.*(1 + x);
                 fclass{1,4} = @(x) 1/4*(-1 + x).*(1 + x).^2;
                 fclass{2,1} = @(x) 3/4*(-1 + x.^2);
                 fclass{2,2} = @(x) -3/4*(-1 + x.^2);
                 fclass{2,3} = @(x) 1/4*(-1 - 2*x + 3*x.^2);
                 fclass{2,4} = @(x) 1/4*(-1 + 2*x + 3*x.^2);
                 fclass{3,1} = @(x) 3*x/2;
                 fclass{3,2} = @(x) -3*x/2;
                 fclass{3,3} = @(x) 1/2*(-1 + 3*x);
                 fclass{3,4} = @(x) 1/2*(1 + 3*x);
                 fclass{4,1} = @(x) 3/2;  
                 fclass{4,2} = @(x) -3/2;  
                 fclass{4,3} = @(x) 3/2;
                 fclass{4,4} = @(x) 3/2;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;  
                 fclass{5,3} = @(x) 0;  fclass{5,4} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 basis.trial.reffun = fclass;
            case 4  % Hermite fourth-order element: a1, da1, a2, a3, da3
                switch femtype.type
                  case 'pgm' 
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                      Pb = zeros(3,mesh.N+1);
                      Pb(1,:) = 1:2:2*mesh.N+1; Pb(2,:) = Pb(1,:);
                      Pb(3,:) = Pb(1,:)+1;
                      Pb = reshape(Pb,1,3*mesh.N+3);
                      Pb = Pb(1:end-1);
                      basis.trial.Pb = Pb0(Pb);
                      Tb = zeros(5,mesh.N);
                      Tb(1,:) = 1:3:3*mesh.N-2; Tb(2,:) = Tb(1,:)+2;
                      Tb(3,:) = Tb(1,:)+3; Tb(4,:) = Tb(1,:)+1;
                      Tb(5,:) = Tb(1,:)+4;
                      basis.trial.Tb = Tb;
                  case 'ipdg' 
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                      Pb = zeros(5,mesh.N);
                      Pb(1,:) = 1:2:2*mesh.N-1;  Pb(2,:) = Pb(1,:);
                      Pb(3,:) = Pb(1,:)+1; Pb(4,:) = Pb(1,:)+2;
                      Pb(5,:) = Pb(1,:)+2;
                      basis.trial.Pb = Pb0(reshape(Pb,1,5*mesh.N));
                      Tb = zeros(5,mesh.N);
                      Tb(1,:) = 1:5:5*mesh.N-4; Tb(2,:) = Tb(1,:)+2;
                      Tb(3,:) = Tb(1,:)+3;  Tb(4,:) = Tb(1,:)+1;
                      Tb(5,:) = Tb(1,:)+4;
                      basis.trial.Tb = Tb;
                end
                 fclass = cell(7,5);
%                  fclass{1,1} = @(x) l1(x)^2*(l1(x)-l2(x))*(l1(x)+5*l2(x));
%                  fclass{1,2} = @(x) 16*l1(x)^2*l2(x)^2;
%                  fclass{1,3} = @(x) -l2(x)^2*(l1(x)-l2(x))*(5*l1(x)+l2(x));
%                  fclass{1,4} = @(x) 2*l1(x)^2*l2(x)*(l1(x)-l2(x));
%                  fclass{1,5} = @(x) 2*l1(x)*l2(x)^2*(l1(x)-l2(x));
                 fclass{1,1} = @(x) -1/4*(-1 + x).^2.*x.*(3 + 2*x);
                 fclass{1,2} = @(x) (-1 + x.^2).^2;
                 fclass{1,3} = @(x) -1/4*x.*(1 + x).^2.*(-3 + 2*x);
                 fclass{1,4} = @(x) -1/4*(-1 + x).^2.*x.*(1 + x);
                 fclass{1,5} = @(x) 1/4*(-1 + x).*x.*(1 + x).^2;
                 fclass{2,1} = @(x) 1/4*(-3 + 8*x).*(1 - x.^2);
                 fclass{2,2} = @(x) 4*x.*(-1 + x.^2);
                 fclass{2,3} = @(x) 1/4*(3 + 8*x).*(1 - x.^2);
                 fclass{2,4} = @(x) 1/4*(-1 + x).*(1 - x - 4*x.^2);
                 fclass{2,5} = @(x) 1/4*(1 + x).*(-1 - x + 4*x.^2);
                 fclass{3,1} = @(x) 2 + 3*x/2 - 6*x.^2;
                 fclass{3,2} = @(x) -4 + 12*x.^2;
                 fclass{3,3} = @(x) 2 - 3*x/2 - 6*x.^2;
                 fclass{3,4} = @(x) 1/2*(1 + 3*x - 6*x.^2);
                 fclass{3,5} = @(x) 1/2*(-1 + 3*x + 6*x.^2);
                 fclass{4,1} = @(x) 3/2 - 12*x;  
                 fclass{4,2} = @(x) 24*x;  
                 fclass{4,3} = @(x) -3/2*(1 + 8*x);
                 fclass{4,4} = @(x) 3/2 - 6*x;
                 fclass{4,5} = @(x) 3/2 + 6*x;
                 fclass{5,1} = @(x) -12;  
                 fclass{5,2} = @(x) 24;  
                 fclass{5,3} = @(x) -12;  
                 fclass{5,4} = @(x) -6; 
                 fclass{5,5} = @(x) 6;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{6,5} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;
                 basis.trial.reffun = fclass;
            case 5  % Hermite fifth-order element: a1, da1, d^2a1, a2, da2, d^2a2
                switch femtype.type
                    case 'pgm' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                        Pb = zeros(3,mesh.N+1);
                        Pb(1,:) = 1:mesh.N+1; Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:);
                        Pb = reshape(Pb,1,3*mesh.N+3);
                        basis.trial.Pb = Pb0(Pb);
                        Tb = zeros(6,mesh.N);
                        Tb(1,:) = 1:3:3*mesh.N-2; Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+1; Tb(4,:) = Tb(1,:)+4;
                        Tb(5,:) = Tb(1,:)+2; Tb(6,:) = Tb(1,:)+5;
                        basis.trial.Tb = Tb;
                    case 'ipdg' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                        Pb = zeros(6,mesh.N);
                        Pb(1,:) = 1:mesh.N;  Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:); Pb(4,:) = Pb(1,:)+1;
                        Pb(5,:) = Pb(1,:)+1; Pb(6,:) = Pb(1,:)+1; 
                        basis.trial.Pb = Pb0(reshape(Pb,1,6*mesh.N));
                        Tb = zeros(6,mesh.N);
                        Tb(1,:) = 1:6:6*mesh.N-5;  Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+1; Tb(4,:) = Tb(1,:)+4;
                        Tb(5,:) = Tb(1,:)+2; Tb(6,:) = Tb(1,:)+5;
                        basis.trial.Tb = Tb;
                end
                 fclass = cell(7,6);
%                  fclass{1,1} = @(x) l1(x)^3*(l1(x)^2+5*l1(x)*l2(x)+10*l2(x)^2);
%                  fclass{1,2} = @(x) l2(x)^3*(10*l1(x)^2+5*l1(x)*l2(x)+l2(x)^2);
%                  fclass{1,3} = @(x) 2*l1(x)^3*l2(x)*(l1(x)+4*l2(x));
%                  fclass{1,4} = @(x) -2*l1(x)*l2(x)^3*(4*l1(x)+l2(x));
%                  fclass{1,5} = @(x) 2*l1(x)^3*l2(x)^2;
%                  fclass{1,6} = @(x) 2*l1(x)^2*l2(x)^3;
                 fclass{1,1} = @(x) -1/16*(-1 + x).^3.*(8 + 9*x + 3*x.^2);
                 fclass{1,2} = @(x) 1/16*(1 + x).^3.*(8 - 9*x + 3*x.^2);
                 fclass{1,3} = @(x) -1/16*(-1 + x).^3.*(1 + x).*(5 + 3*x);
                 fclass{1,4} = @(x) -1/16*(-1 + x).*(1 + x).^3.*(-5 + 3*x);
                 fclass{1,5} = @(x) -1/16*(-1 + x).^3.*(1 + x).^2;
                 fclass{1,6} = @(x) 1/16*(-1 + x).^2.*(1 + x).^3;
                 fclass{2,1} = @(x) -15/16*(-1 + x.^2).^2;
                 fclass{2,2} = @(x) 15/16*(-1 + x.^2).^2;
                 fclass{2,3} = @(x) -1/16*(-1 + x).^2.*(1 + 3*x).*(7 + 5*x);
                 fclass{2,4} = @(x) -1/16*(1 + x).^2.*(-1 + 3*x).*(-7 + 5*x);
                 fclass{2,5} = @(x) -1/16*(-1 + x).^2.*(1 + x).*(1 + 5*x);
                 fclass{2,6} = @(x) 1/16*(1 + x).^2.*(-1 + x).*(-1 + 5*x);
                 fclass{3,1} = @(x) -15/4*x.*(-1 + x.^2);
                 fclass{3,2} = @(x) 15/4*x.*(-1 + x.^2);
                 fclass{3,3} = @(x) -3/4*(-1 + 5*x).*(- 1 + x.^2);
                 fclass{3,4} = @(x) -3/4*(1 + 5*x).*(- 1 + x.^2);
                 fclass{3,5} = @(x) 1/4*(-1 + 3*x + 3*x.^2 - 5*x.^3);
                 fclass{3,6} = @(x) 1/4*(-1 - 3*x + 3*x.^2 + 5*x.^3);
                 fclass{4,1} = @(x) -15/4*(-1 + 3*x.^2); 
                 fclass{4,2} = @(x) 15/4*(-1 + 3*x.^2);  
                 fclass{4,3} = @(x) -3/4*(-5 - 2*x + 15*x.^2);
                 fclass{4,4} = @(x) -3/4*(-5 + 2*x + 15*x.^2);
                 fclass{4,5} = @(x) -3/4*(-1 - 2*x + 5*x.^2);
                 fclass{4,6} = @(x) 3/4*(-1 + 2*x + 5*x.^2); 
                 fclass{5,1} = @(x) -45*x/2;  
                 fclass{5,2} = @(x) 45*x/2;  
                 fclass{5,3} = @(x) 3/2 - 45*x/2;  
                 fclass{5,4} = @(x) -3/2*(1 + 15*x); 
                 fclass{5,5} = @(x) -3/2*(-1 + 5*x);
                 fclass{5,6} = @(x) 3/2*(1 + 5*x);
                 fclass{6,1} = @(x) -45/2;  
                 fclass{6,2} = @(x) 45/2;  
                 fclass{6,3} = @(x) -45/2;  
                 fclass{6,4} = @(x) -45/2;
                 fclass{6,5} = @(x) -15/2;  
                 fclass{6,6} = @(x) 15/2;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;  fclass{7,6} = @(x) 0;
                 basis.trial.reffun = fclass;
            case 6  % Hermite sixth-order element: a1, da1, d^2a1, a2, a3, da3, d^2a3
                switch femtype.type
                    case 'pgm' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                        Pb = zeros(4,mesh.N+1);
                        Pb(1,:) = 1:2:2*mesh.N+1; Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:); Pb(4,:) = Pb(1,:)+1;
                        Pb = Pb(1:end-1);
                        basis.trial.Pb = Pb0(reshape(Pb,1,4*mesh.N+3));
                        Tb = zeros(7,mesh.N);
                        Tb(1,:) = 1:4:4*mesh.N-3; Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+4; Tb(4,:) = Tb(1,:)+1;
                        Tb(5,:) = Tb(1,:)+5; Tb(6,:) = Tb(1,:)+2;
                        Tb(7,:) = Tb(1,:)+6;
                        basis.trial.Tb = Tb;
                    case 'ipdg' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                        Pb = zeros(7,mesh.N);
                        Pb(1,:) = 1:2:2*mesh.N-1;  Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:); Pb(4,:) = Pb(1,:)+1;
                        Pb(5,:) = Pb(1,:)+2; Pb(6,:) = Pb(1,:)+2; 
                        Pb(7,:) = Pb(1,:)+2;
                        basis.trial.Pb = Pb0(reshape(Pb,1,7*mesh.N));
                        Tb = zeros(7,mesh.N);
                        Tb(1,:) = 1:7:7*mesh.N-6;  Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+4; Tb(4,:) = Tb(1,:)+1;
                        Tb(5,:) = Tb(1,:)+5; Tb(6,:) = Tb(1,:)+2;
                        Tb(7,:) = Tb(1,:)+6;
                        basis.trial.Tb = Tb;
                end
                 fclass = cell(7,7);
%                  fclass{1,1} = @(x) l1(x)^3*(l1(x)-l2(x))*(l1(x)^2+7*l1(x)*l2(x)+22*l2(x)^2);
%                  fclass{1,2} = @(x) 64*l1(x)^3*l2(x)^3;
%                  fclass{1,3} = @(x) -l2(x)^3*(l1(x)-l2(x))*(22*l1(x)^2+7*l1(x)*l2(x)+l2(x)^2);
%                  fclass{1,4} = @(x) l1(x)^3*l2(x)*(l1(x)-l2(x))*(2*l1(x)+12*l2(x));
%                  fclass{1,5} = @(x) l1(x)*l2(x)^3*(l1(x)-l2(x))*(12*l1(x)+2*l2(x));
%                  fclass{1,6} = @(x) 2*l1(x)^3*l2(x)^2*(l1(x)-l2(x));
%                  fclass{1,7} = @(x) -2*l1(x)^2*l2(x)^3*(l1(x)-l2(x));
                 fclass{1,1} = @(x) 1/16*(-1 + x).^3.*x.*(15 + 21*x + 8*x.^2);
                 fclass{1,2} = @(x) -(-1 + x.^2).^3;
                 fclass{1,3} = @(x) 1/16*x.*(1 + x).^3.*(15 - 21*x + 8*x.^2);
                 fclass{1,4} = @(x) 1/16*(-1 + x).^3.*x.*(1 + x).*(7 + 5*x);
                 fclass{1,5} = @(x) -1/16*(-1 + x).*x.*(1 + x).^3.*(-7 + 5*x);
                 fclass{1,6} = @(x) 1/16*(-1 + x).^3.*x.*(1 + x).^2;
                 fclass{1,7} = @(x) 1/16*(-1 + x).^2.*x.*(1 + x).^3;
                 fclass{2,1} = @(x) 3/16*(-5 + 16*x).*(-1 + x.^2).^2;
                 fclass{2,2} = @(x) -6*x.*(-1 + x.^2).^2;
                 fclass{2,3} = @(x) 3/16*(5 + 16*x).*(-1 + x.^2).^2;
                 fclass{2,4} = @(x) 1/16*(-1 + x).^2.*(-7 + 4*x + 45*x.^2 + 30*x.^3);
                 fclass{2,5} = @(x) -1/16*(1 + x).^2.*(7 + 4*x - 45*x.^2 + 30*x.^3);
                 fclass{2,6} = @(x) 1/16*(-1 + x).^2.*(-1 + 7*x.^2 + 6*x.^3);
                 fclass{2,7} = @(x) 1/16*(1 + x).^2.*(1 - 7*x.^2 + 6*x.^3);
                 fclass{3,1} = @(x) 3/4*(4 + 5*x - 24*x.^2 - 5*x.^3 + 20*x.^4);
                 fclass{3,2} = @(x) -6*(1 - 6*x.^2 + 5*x.^4);
                 fclass{3,3} = @(x) 3/4*(4 - 5*x - 24*x.^2 + 5*x.^3 + 20*x.^4);
                 fclass{3,4} = @(x) 3/8*(3 + 10*x - 28*x.^2 - 10*x.^3 + 25*x.^4);
                 fclass{3,5} = @(x) -3/8*(3 - 10*x - 28*x.^2 + 10*x.^3 + 25*x.^4);
                 fclass{3,6} = @(x) 1/8*(1 + 6*x - 12*x.^2 - 10*x.^3 + 15*x.^4);
                 fclass{3,7} = @(x) 1/8*(1 - 6*x - 12*x.^2 + 10*x.^3 + 15*x.^4);
                 fclass{4,1} = @(x) 3/4*(5 - 48*x - 15*x.^2 + 80*x.^3);  
                 fclass{4,2} = @(x) 24*x.*(3 - 5*x.^2);  
                 fclass{4,3} = @(x) 3/4*(-5 - 48*x + 15*x.^2 + 80*x.^3);
                 fclass{4,4} = @(x) 3/4*(5 - 28*x - 15*x.^2 + 50*x.^3);
                 fclass{4,5} = @(x) -3/4*(-5 - 28*x + 15*x.^2 + 50*x.^3);
                 fclass{4,6} = @(x) 3/4*(1 - 4*x - 5*x.^2 + 10*x.^3);
                 fclass{4,7} = @(x) 3/4*(-1 - 4*x + 5*x.^2 + 10*x.^3);
                 fclass{5,1} = @(x) 9/2*(-8 - 5*x + 40*x.^2);  
                 fclass{5,2} = @(x) 72 - 360*x.^2;  
                 fclass{5,3} = @(x) 9/2*(-8 + 5*x + 40*x.^2);  
                 fclass{5,4} = @(x) 3/2*(-14 - 15*x + 75*x.^2); 
                 fclass{5,5} = @(x) -3/2*(-14 + 15*x + 75*x.^2);
                 fclass{5,6} = @(x) 3/2*(-2 - 5*x + 15*x.^2);
                 fclass{5,7} = @(x) 3/2*(-2 + 5*x + 15*x.^2);
                 fclass{6,1} = @(x) -45/2 + 360*x;  
                 fclass{6,2} = @(x) -720*x;  
                 fclass{6,3} = @(x) 45/2 + 360*x;  
                 fclass{6,4} = @(x) -45/2 + 225*x;
                 fclass{6,5} = @(x) -45/2*(1 + 10*x);  
                 fclass{6,6} = @(x) -15/2 + 45*x;
                 fclass{6,7} = @(x) 15/2 + 45*x;
                 fclass{7,1} = @(x) 360;  
                 fclass{7,2} = @(x) -720;  
                 fclass{7,3} = @(x) 360;  
                 fclass{7,4} = @(x) 225;
                 fclass{7,5} = @(x) -450/2;  
                 fclass{7,6} = @(x) 45;
                 fclass{7,7} = @(x) 45;
                 basis.trial.reffun = fclass;   
        end
    basis.trial.Nb = size(basis.trial.Pb,2);
    basis.trial.Nlb = size(basis.trial.Tb,1);
end

%% generate test basis function as well as Pb and Tb matrix, 
% the basis function is defined on the reference element (0,1)
if basis.test.type ~= 'l' && basis.test.type ~= 'h'
    warning('The basis.test.type must be either l or h, already reset to l');
    basis.test.type = 'l';
end
if basis.test.type == 'l'
    ord = basis.test.order;
    if floor(ord)~=ord || ord<1 || ord>6
        warning(['The basis.test.order must be a positive integer less than 7,',...
                'already reset to 1']);
        basis.test.order = 1;
    end
    switch basis.test.order
          case 1 % Lagrange linear basis function
              switch femtype.type
                  case 'pgm'
                 basis.test.Pb = mesh.P;
                 basis.test.Tb = mesh.T;
                  case 'ipdg'
                 Pb = mesh.P;
                 Pb2 = Pb(:,2:end-1);
                 basis.test.Pb = [Pb(:,1),reshape([Pb2;Pb2],1,2*mesh.N-2),Pb(:,end)];
                 basis.test.Tb = reshape(1:2*mesh.N,2,mesh.N);     
              end
                 fclass = cell(7,2);
%                  fclass{1,1} = @(x) l1(x);
%                  fclass{1,2} = @(x) l2(x);
                 fclass{1,1} = @(x) 1/2*(-x+1);
                 fclass{1,2} = @(x) 1/2*(x+1);
                 fclass{2,1} = @(x) -1/2;
                 fclass{2,2} = @(x) 1/2;
                 fclass{3,1} = @(x) 0;  fclass{3,2} = @(x) 0;
                 fclass{4,1} = @(x) 0;  fclass{4,2} = @(x) 0;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;
                 basis.test.reffun = fclass;
          case 2 % Lagrange quadratic basis function
              switch femtype.type
                  case 'pgm'
                      basis.test.Pb = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                      Tb = zeros(3,mesh.N);
                      Tb(1,:) = 1:2:2*mesh.N-1; Tb(2,:) = 2:2:2*mesh.N;
                      Tb(3,:) = 3:2:2*mesh.N+1;
                      basis.test.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                      Pb = zeros(3,mesh.N);
                      Pb(1,:) = 1:2:2*mesh.N-1; Pb(2,:) = 2:2:2*mesh.N;
                      Pb(3,:) = 3:2:2*mesh.N+1;
                      Pb = reshape(Pb,1,3*mesh.N);
                      basis.test.Pb = Pb0(Pb);
                      basis.test.Tb = reshape(1:3*mesh.N,3,mesh.N);
              end
                 fclass = cell(7,3);
%                  fclass{1,1} = @(x) l1(x)*(l1(x)-l2(x));
%                  fclass{1,2} = @(x) 4*l1(x)*l2(x);
%                  fclass{1,3} = @(x) l2(x)*(l2(x)-l1(x));
                 fclass{1,1} = @(x) 1/2*(-1 + x).*x;
                 fclass{1,2} = @(x) 1-x.^2;
                 fclass{1,3} = @(x) 1/2*x.*(1 + x);
                 fclass{2,1} = @(x) -1/2 + x;
                 fclass{2,2} = @(x) -2*x;
                 fclass{2,3} = @(x) 1/2 + x;
                 fclass{3,1} = @(x) 1;
                 fclass{3,2} = @(x) -2;
                 fclass{3,3} = @(x) 1;
                 fclass{4,1} = @(x) 0;  fclass{4,2} = @(x) 0;  fclass{4,3} = @(x) 0;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;  fclass{5,3} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  fclass{6,3} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  fclass{7,3} = @(x) 0;
                 basis.test.reffun = fclass;
          case 3 % Lagrange cubic basis function
              switch femtype.type
                  case 'pgm'
                      basis.test.Pb = mesh.left:(mesh.right-mesh.left)/(3*mesh.N):mesh.right;
                      Tb = zeros(4,mesh.N);
                      Tb(1,:) = 1:3:3*mesh.N-2; Tb(2,:) = 2:3:3*mesh.N-1; 
                      Tb(3,:) = 3:3:3*mesh.N; Tb(4,:) = 4:3:3*mesh.N+1;
                      basis.test.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(3*mesh.N):mesh.right;
                      Pb = zeros(4,mesh.N);
                      Pb(1,:) = 1:3:3*mesh.N-2; Pb(2,:) = 2:3:3*mesh.N-1; 
                      Pb(3,:) = 3:3:3*mesh.N; Pb(4,:) = 4:3:3*mesh.N+1;
                      Pb = reshape(Pb,1,4*mesh.N);
                      basis.test.Pb = Pb0(Pb);
                      basis.test.Tb = reshape(1:4*mesh.N,4,mesh.N);
               end
                 fclass = cell(7,4);
%                  fclass{1,1} = @(x) 1/2*l1(x)*(l1(x)-2*l2(x))*(2*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 9/2*l1(x)*l2(x)*(2*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -9/2*l1(x)*l2(x)*(l1(x)-2*l2(x));
%                  fclass{1,4} = @(x) 1/2*l2(x)*(l1(x)-2*l2(x))*(2*l1(x)-l2(x));
                 fclass{1,1} = @(x) -1/16*(-1 + x).*(-1 + 3*x).*(1 + 3*x);
                 fclass{1,2} = @(x) 9/16*(-1 + x).*(1 + x).*(-1 + 3*x);
                 fclass{1,3} = @(x) -9/16*(-1 + x).*(1 + x).*(1 + 3*x);
                 fclass{1,4} = @(x) 1/16*(1 + x).*(-1 + 3*x).*(1 + 3*x);
                 fclass{2,1} = @(x) 1/16*(1 + 18*x - 27*x.^2);
                 fclass{2,2} = @(x) 9/16*(-3 - 2*x + 9*x.^2);
                 fclass{2,3} = @(x) -9/16*(-3 + 2*x + 9*x.^2);
                 fclass{2,4} = @(x) 1/16*(-1 + 18*x + 27*x.^2);
                 fclass{3,1} = @(x) -9/8*(-1 + 3*x);
                 fclass{3,2} = @(x) 9/8*(-1 + 9*x);
                 fclass{3,3} = @(x) -9/8*(1 + 9*x);
                 fclass{3,4} = @(x) 9/8*(1 + 3*x);
                 fclass{4,1} = @(x) -27/8;  
                 fclass{4,2} = @(x) 81/8;  
                 fclass{4,3} = @(x) -81/8;
                 fclass{4,4} = @(x) 27/8;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;  
                 fclass{5,3} = @(x) 0;  fclass{5,4} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 basis.test.reffun = fclass;
          case 4 % Lagrange fourth-order basis function
              switch femtype.type
                  case 'pgm'
                      basis.test.Pb = mesh.left:(mesh.right-mesh.left)/(4*mesh.N):mesh.right;
                      Tb = zeros(5,mesh.N);
                      Tb(1,:) = 1:4:4*mesh.N-3; Tb(2,:) = 2:4:4*mesh.N-2; 
                      Tb(3,:) = 3:4:4*mesh.N-1; Tb(4,:) = 4:4:4*mesh.N;
                      Tb(5,:) = 5:4:4*mesh.N+1;
                      basis.test.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(4*mesh.N):mesh.right;
                      Pb = zeros(5,mesh.N);
                      Pb(1,:) = 1:4:4*mesh.N-3; Pb(2,:) = 2:4:4*mesh.N-2; 
                      Pb(3,:) = 3:4:4*mesh.N-1; Pb(4,:) = 4:4:4*mesh.N;
                      Pb(5,:) = 5:4:4*mesh.N+1;
                      Pb = reshape(Pb,1,5*mesh.N);
                      basis.test.Pb = Pb0(Pb);
                      basis.test.Tb = reshape(1:5*mesh.N,5,mesh.N);
              end
                 fclass = cell(7,5);
%                  fclass{1,1} = @(x) 1/3*l1(x)*(l1(x)-l2(x))*(l1(x)-3*l2(x))*(3*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 16/3*l1(x)*l2(x)*(l1(x)-l2(x))*(3*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -4*l1(x)*l2(x)*(l1(x)-3*l2(x))*(3*l1(x)-l2(x));
%                  fclass{1,4} = @(x) 16/3*l1(x)*l2(x)*(l1(x)-l2(x))*(l1(x)-3*l2(x));
%                  fclass{1,5} = @(x) -1/3*l2(x)*(l1(x)-l2(x))*(l1(x)-3*l2(x))*(3*l1(x)-l2(x));
                 fclass{1,1} = @(x) 1/6*(-1 + x).*x.*(-1 + 2*x).*(1 + 2*x);
                 fclass{1,2} = @(x) -4/3*(-1 + x).*x.*(1 + x).*(-1 + 2*x);
                 fclass{1,3} = @(x) 1 - 5*x.^2 + 4*x.^4;
                 fclass{1,4} = @(x) -4/3*(-1 + x).*x.*(1 + x).*(1 + 2*x);
                 fclass{1,5} = @(x) 1/6*x.*(1 + x).*(-1 + 2*x).*(1 + 2*x);
                 fclass{2,1} = @(x) 1/6*(1 - 2*x - 12*x.^2 + 16*x.^3);
                 fclass{2,2} = @(x) -4/3*(1 - 4*x - 3*x.^2 + 8*x.^3);
                 fclass{2,3} = @(x) 2*x.*(-5 + 8*x.^2);
                 fclass{2,4} = @(x) -4/3*(-1 - 4*x + 3*x.^2 + 8*x.^3);
                 fclass{2,5} = @(x) 1/6*(-1 - 2*x + 12*x.^2 + 16*x.^3);
                 fclass{3,1} = @(x) -1/3 - 4*x + 8*x.^2;
                 fclass{3,2} = @(x) 16/3 + 8*x - 32*x.^2;
                 fclass{3,3} = @(x) -10 + 48*x.^2;
                 fclass{3,4} = @(x) 16/3 - 8*x - 32*x.^2;
                 fclass{3,5} = @(x) -1/3 + 4*x + 8*x.^2;
                 fclass{4,1} = @(x) -4 + 16*x;  
                 fclass{4,2} = @(x) 8 - 64*x;  
                 fclass{4,3} = @(x) 96*x;
                 fclass{4,4} = @(x) -8*(1 + 8*x);
                 fclass{4,5} = @(x) 4 + 16*x;
                 fclass{5,1} = @(x) 16;  
                 fclass{5,2} = @(x) -64;  
                 fclass{5,3} = @(x) 96;  
                 fclass{5,4} = @(x) -64; 
                 fclass{5,5} = @(x) 16;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{6,5} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;
                 basis.test.reffun = fclass;
          case 5 % Lagrange fifth-order basis function
              switch femtype.type
                  case 'pgm'
                      basis.test.Pb = mesh.left:(mesh.right-mesh.left)/(5*mesh.N):mesh.right;
                      Tb = zeros(6,mesh.N);
                      Tb(1,:) = 1:5:5*mesh.N-4; Tb(2,:) = 2:5:5*mesh.N-3; 
                      Tb(3,:) = 3:5:5*mesh.N-2; Tb(4,:) = 4:5:5*mesh.N-1;
                      Tb(5,:) = 5:5:5*mesh.N; Tb(6,:) = 6:5:5*mesh.N+1;
                      basis.test.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(5*mesh.N):mesh.right;
                      Pb = zeros(6,mesh.N);
                      Pb(1,:) = 1:5:5*mesh.N-4; Pb(2,:) = 2:5:5*mesh.N-3; 
                      Pb(3,:) = 3:5:5*mesh.N-2; Pb(4,:) = 4:5:5*mesh.N-1;
                      Pb(5,:) = 5:5:5*mesh.N; Pb(6,:) = 6:5:5*mesh.N+1;
                      Pb = reshape(Pb,1,6*mesh.N);
                      basis.test.Pb = Pb0(Pb);
                      basis.test.Tb = reshape(1:6*mesh.N,6,mesh.N);
              end
                 fclass = cell(7,6);
%                  fclass{1,1} = @(x) 1/24*l1(x)*(l1(x)-4*l2(x))*(2*l1(x)-3*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 25/24*l1(x)*l2(x)*(2*l1(x)-3*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -25/12*l1(x)*l2(x)*(l1(x)-4*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,4} = @(x) 25/12*l1(x)*l2(x)*(l1(x)-4*l2(x))...
%                                     *(2*l1(x)-3*l2(x))*(4*l1(x)-l2(x));
%                  fclass{1,5} = @(x) -25/24*l1(x)*l2(x)*(l1(x)-4*l2(x))...
%                                     *(2*l1(x)-3*l2(x))*(3*l1(x)-2*l2(x));
%                  fclass{1,6} = @(x) 1/24*l2(x)*(l1(x)-4*l2(x))*(2*l1(x)-3*l2(x))...
%                                     *(3*l1(x)-2*l2(x))*(4*l1(x)-l2(x));
                 fclass{1,1} = @(x) -1/768*(-1 + x).*(-3 + 5*x).*(-1 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{1,2} = @(x) 25/768*(-1 + x).*(1 + x).*(-3 + 5*x).*(-1 + 5*x).*(1 + 5*x);
                 fclass{1,3} = @(x) -25/384*(-1 + x).*(1 + x).*(-3 + 5*x).*(-1 + 5*x).*(3 + 5*x);
                 fclass{1,4} = @(x) 25/384*(-1 + x).*(1 + x).*(-3 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{1,5} = @(x) -25/768*(-1 + x).*(1 + x).*(-1 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{1,6} = @(x) 1/768*(1 + x).*(-3 + 5*x).*(-1 + 5*x).*(1 + 5*x).*(3 + 5*x);
                 fclass{2,1} = @(x) 1/768*(-9 - 500*x + 750*x.^2 + 2500*x.^3 - 3125*x.^4);
                 fclass{2,2} = @(x) 25/768*(5 + 156*x - 390*x.^2 - 300*x.^3 + 625*x.^4);
                 fclass{2,3} = @(x) -25/384*(45 + 68*x - 510*x.^2 - 100*x.^3 + 625*x.^4);
                 fclass{2,4} = @(x) 25/384*(45 - 68*x - 510*x.^2 + 100*x.^3 + 625*x.^4);
                 fclass{2,5} = @(x) -25/768*(5 - 156*x - 390*x.^2 + 300*x.^3 + 625*x.^4);
                 fclass{2,6} = @(x) 1/768*(9 - 500*x - 750*x.^2 + 2500*x.^3 + 3125*x.^4);
                 fclass{3,1} = @(x) -125/192*(1 - 3*x - 15*x.^2 + 25*x.^3);
                 fclass{3,2} = @(x) 25/192*(39 - 195*x - 225*x.^2 + 625*x.^3);
                 fclass{3,3} = @(x) -25/96*(17 - 255*x - 75*x.^2 + 625*x.^3);
                 fclass{3,4} = @(x) 25/96*(-17 - 255*x + 75*x.^2 + 625*x.^3);
                 fclass{3,5} = @(x) -25/192*(-39 - 195*x + 225*x.^2 + 625*x.^3);
                 fclass{3,6} = @(x) 125/192*(-1 - 3*x + 15*x.^2 + 25*x.^3);
                 fclass{4,1} = @(x) -125/64*(-1 - 10*x + 25*x.^2);  
                 fclass{4,2} = @(x) 125/64*(-13 - 30*x + 125*x.^2);  
                 fclass{4,3} = @(x) -125/32*(-17 - 10*x + 125*x.^2);
                 fclass{4,4} = @(x) 125/32*(-17 + 10*x + 125*x.^2);
                 fclass{4,5} = @(x) -125/64*(-13 + 30*x + 125*x.^2);
                 fclass{4,6} = @(x) 125/64*(-1 + 10*x + 25*x.^2);
                 fclass{5,1} = @(x) -625/32*(-1 + 5*x);  
                 fclass{5,2} = @(x) 625/32*(-3 + 25*x);  
                 fclass{5,3} = @(x) -625/16*(-1 + 25*x);  
                 fclass{5,4} = @(x) 625/16*(1 + 25*x); 
                 fclass{5,5} = @(x) -625/32*(3 + 25*x);
                 fclass{5,6} = @(x) 625/32*(1 + 5*x);
                 fclass{6,1} = @(x) -3125/32;  
                 fclass{6,2} = @(x) 15625/32;  
                 fclass{6,3} = @(x) -15625/16;  
                 fclass{6,4} = @(x) 15625/16;
                 fclass{6,5} = @(x) -15625/32;  
                 fclass{6,6} = @(x) 3125/32;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;  fclass{7,6} = @(x) 0;
                 basis.test.reffun = fclass;
          case 6 % Lagrange sixth-order basis function
              switch femtype.type
                  case 'pgm'
                      basis.test.Pb = mesh.left:(mesh.right-mesh.left)/(6*mesh.N):mesh.right;
                      Tb = zeros(7,mesh.N);
                      Tb(1,:) = 1:6:6*mesh.N-5; Tb(2,:) = 2:6:6*mesh.N-4; 
                      Tb(3,:) = 3:6:6*mesh.N-3; Tb(4,:) = 4:6:6*mesh.N-2;
                      Tb(5,:) = 5:6:6*mesh.N-1; Tb(6,:) = 6:6:6*mesh.N;
                      Tb(7,:) = 7:6:6*mesh.N+1;
                      basis.test.Tb = Tb;
                  case 'ipdg'
                      Pb0 = mesh.left:(mesh.right-mesh.left)/(6*mesh.N):mesh.right;
                      Pb = zeros(7,mesh.N);
                      Pb(1,:) = 1:6:6*mesh.N-5; Pb(2,:) = 2:6:6*mesh.N-4; 
                      Pb(3,:) = 3:6:6*mesh.N-3; Pb(4,:) = 4:6:6*mesh.N-2;
                      Pb(5,:) = 5:6:6*mesh.N-1; Pb(6,:) = 6:6:6*mesh.N;
                      Pb(7,:) = 7:6:6*mesh.N+1;
                      Pb = reshape(Pb,1,7*mesh.N);
                      basis.test.Pb = Pb0(Pb);
                      basis.test.Tb = reshape(1:7*mesh.N,7,mesh.N);
              end
                 fclass = cell(7,7);
%                  fclass{1,1} = @(x) 1/10*l1(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,2} = @(x) 18/5*l1(x)*l2(x)*(l1(x)-2*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,3} = @(x) -9/2*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,4} = @(x) 4*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,5} = @(x) -9/2*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))...
%                                     *(l1(x)-l2(x))*(5*l1(x)-l2(x));
%                  fclass{1,6} = @(x) 18/5*l1(x)*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))...
%                                     *(l1(x)-l2(x))*(2*l1(x)-l2(x));
%                  fclass{1,7} = @(x) -1/10*l2(x)*(l1(x)-5*l2(x))*(l1(x)-2*l2(x))*(l1(x)-l2(x))...
%                                     *(2*l1(x)-l2(x))*(5*l1(x)-l2(x));
                 fclass{1,1} = @(x) 1/80*(-1 + x).*x.*(-2 + 3*x).*(-1 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{1,2} = @(x) -9/40*(-1 + x).*x.*(1 + x).*(-2 + 3*x).*(-1 + 3*x).*(1 + 3*x);
                 fclass{1,3} = @(x) 9/16*(-1 + x).*x.*(1 + x).*(-2 + 3*x).*(-1 + 3*x).*(2 + 3*x);
                 fclass{1,4} = @(x) 1/4*(4 - 49*x.^2 + 126*x.^4 - 81*x.^6);
                 fclass{1,5} = @(x) 9/16*(-1 + x).*x.*(1 + x).*(-2 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{1,6} = @(x) -9/40*(-1 + x).*x.*(1 + x).*(-1 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{1,7} = @(x) 1/80*x.*(1 + x).*(-2 + 3*x).*(-1 + 3*x).*(1 + 3*x).*(2 + 3*x);
                 fclass{2,1} = @(x) 1/80*(-4 + 8*x + 135*x.^2 - 180*x.^3 - 405*x.^4 + 486*x.^5);
                 fclass{2,2} = @(x) -9/20*(-1 + 3*x + 30*x.^2 - 60*x.^3 - 45*x.^4 + 81*x.^5);
                 fclass{2,3} = @(x) 9/16*(-4 + 24*x + 39*x.^2 - 156*x.^3 - 45*x.^4 + 162*x.^5);
                 fclass{2,4} = @(x) 1/2*(-49*x + 252*x.^3 - 243*x.^5);
                 fclass{2,5} = @(x) 9/16*(4 + 24*x - 39*x.^2 - 156*x.^3 + 45*x.^4 + 162*x.^5);
                 fclass{2,6} = @(x) -9/20*(1 + 3*x - 30*x.^2 - 60*x.^3 + 45*x.^4 + 81*x.^5);
                 fclass{2,7} = @(x) 1/80*(4 + 8*x - 135*x.^2 - 180*x.^3 + 405*x.^4 + 486*x.^5);
                 fclass{3,1} = @(x) 1/40*(4 + 135*x - 270*x.^2 - 810*x.^3 + 1215*x.^4);
                 fclass{3,2} = @(x) -27/20*(1 + 20*x - 60*x.^2 - 60*x.^3 + 135*x.^4);
                 fclass{3,3} = @(x) 27/8*(4 + 13*x - 78*x.^2 - 30*x.^3 + 135*x.^4);
                 fclass{3,4} = @(x) 1/2*(-49 + 756*x.^2 - 1215*x.^4);
                 fclass{3,5} = @(x) 27/8*(4 - 13*x - 78*x.^2 + 30*x.^3 + 135*x.^4);
                 fclass{3,6} = @(x) -27/20*(1 - 20*x - 60*x.^2 + 60*x.^3 + 135*x.^4);
                 fclass{3,7} = @(x) 1/40*(4 - 135*x - 270*x.^2 + 810*x.^3 + 1215*x.^4);
                 fclass{4,1} = @(x) 27/8*(1 - 4*x - 18*x.^2 + 36*x.^3);  
                 fclass{4,2} = @(x) -27*(1 - 6*x - 9*x.^2 + 27*x.^3);  
                 fclass{4,3} = @(x) 27/8*(13 - 156*x - 90*x.^2 + 540*x.^3);
                 fclass{4,4} = @(x) 54*x.*(14 - 45*x.^2);
                 fclass{4,5} = @(x) 27/8*(-13 - 156*x + 90*x.^2 + 540*x.^3);
                 fclass{4,6} = @(x) -27*(-1 - 6*x + 9*x.^2 + 27.*x.^3);
                 fclass{4,7} = @(x) 27/8*(-1 - 4*x + 18*x.^2 + 36*x.^3);
                 fclass{5,1} = @(x) 27/2*(-1 - 9*x + 27*x.^2);  
                 fclass{5,2} = @(x) -81*(-2 - 6*x + 27*x.^2);  
                 fclass{5,3} = @(x) 81/2*(-13 - 15*x + 135*x.^2);  
                 fclass{5,4} = @(x) 756 - 7290*x.^2; 
                 fclass{5,5} = @(x) 81/2*(-13 + 15*x + 135*x.^2);
                 fclass{5,6} = @(x) -81*(-2 + 6*x + 27*x.^2);
                 fclass{5,7} = @(x) 27/2*(-1 + 9*x + 27*x.^2);
                 fclass{6,1} = @(x) 243/2*(-1 + 6*x);  
                 fclass{6,2} = @(x) -486*(-1 + 9*x);  
                 fclass{6,3} = @(x) 1215/2*(-1 + 18*x);  
                 fclass{6,4} = @(x) -14580*x;
                 fclass{6,5} = @(x) 1215/2*(1 + 18*x);  
                 fclass{6,6} = @(x) -486*(1 + 9*x);
                 fclass{6,7} = @(x) 243/2*(1 + 6*x);
                 fclass{7,1} = @(x) 729;  
                 fclass{7,2} = @(x) -4374;  
                 fclass{7,3} = @(x) 10935;  
                 fclass{7,4} = @(x) -14580;
                 fclass{7,5} = @(x) 10935;  
                 fclass{7,6} = @(x) -4374;
                 fclass{7,7} = @(x) 729;
                 basis.test.reffun = fclass;
    end
    basis.test.Nb = size(basis.test.Pb,2);
    basis.test.Nlb = size(basis.test.Tb,1);
elseif basis.test.type == 'h'
    ord = basis.test.order;
    if ord<3 || ord>6
        warning(['The basis.test.order must be an interger more than 2 less than 7 ,',...
                'already reset to 3']);
        basis.trial.order = 3;
    end
        switch basis.test.order 
            case 3  % Hermite cubic element: a1, da1, a2, da2
               switch femtype.type
                  case 'pgm' 
                      Pb = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                      basis.test.Pb = reshape([Pb; Pb],1,2*mesh.N+2);
                      Tb = zeros(4,mesh.N);
                      Tb(1,:) = 1:2:2*mesh.N-1; Tb(2,:) = Tb(1,:)+2; 
                      Tb(3,:) = Tb(1,:)+1; Tb(4,:) = Tb(1,:)+3; 
                      basis.test.Tb = Tb;
                  case 'ipdg' 
                      Pb0 = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                      Tb = reshape(1:4*mesh.N,4,mesh.N);
                      t = Tb(2,:); Tb(2,:) = Tb(3,:); Tb(3,:) = t;
                      basis.test.Tb = Tb;
                      Pb = zeros(4,mesh.N);
                      Pb(1,:) = 1:mesh.N;  Pb(2,:) = Pb(1,:);
                      Pb(3,:) = 2:mesh.N+1; Pb(4,:) = Pb(3,:);
                      basis.test.Pb = Pb0(reshape(Pb,1,4*mesh.N));
               end
                 fclass = cell(7,4);
%                  fclass{1,1} = @(x) l1(x)^2*(l1(x)+3*l2(x));
%                  fclass{1,2} = @(x) l2(x)^2*(3*l1(x)+l2(x));
%                  fclass{1,3} = @(x) 2*l1(x)^2*l2(x);
%                  fclass{1,4} = @(x) -2*l1(x)*l2(x)^2;
                 fclass{1,1} = @(x) 1/4*(-1 + x).^2.*(2 + x);
                 fclass{1,2} = @(x) -1/4*(-2 + x).*(1 + x).^2;
                 fclass{1,3} = @(x) 1/4*(-1 + x).^2.*(1 + x);
                 fclass{1,4} = @(x) 1/4*(-1 + x).*(1 + x).^2;
                 fclass{2,1} = @(x) 3/4*(-1 + x.^2);
                 fclass{2,2} = @(x) -3/4*(-1 + x.^2);
                 fclass{2,3} = @(x) 1/4*(-1 - 2*x + 3*x.^2);
                 fclass{2,4} = @(x) 1/4*(-1 + 2*x + 3*x.^2);
                 fclass{3,1} = @(x) 3*x/2;
                 fclass{3,2} = @(x) -3*x/2;
                 fclass{3,3} = @(x) 1/2*(-1 + 3*x);
                 fclass{3,4} = @(x) 1/2*(1 + 3*x);
                 fclass{4,1} = @(x) 3/2;  
                 fclass{4,2} = @(x) -3/2;  
                 fclass{4,3} = @(x) 3/2;
                 fclass{4,4} = @(x) 3/2;
                 fclass{5,1} = @(x) 0;  fclass{5,2} = @(x) 0;  
                 fclass{5,3} = @(x) 0;  fclass{5,4} = @(x) 0;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 basis.test.reffun = fclass;
            case 4  % Hermite fourth-order element: a1, da1, a2, a3, da3
                switch femtype.type
                    case 'pgm' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                        Pb = zeros(3,mesh.N+1);
                        Pb(1,:) = 1:2:2*mesh.N+1; Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:)+1;
                        Pb = reshape(Pb,1,3*mesh.N+3);
                        Pb = Pb(1:end-1);
                        basis.test.Pb = Pb0(Pb);
                        Tb = zeros(5,mesh.N);
                        Tb(1,:) = 1:3:3*mesh.N-2; Tb(2,:) = Tb(1,:)+2;
                        Tb(3,:) = Tb(1,:)+3; Tb(4,:) = Tb(1,:)+1;
                        Tb(5,:) = Tb(1,:)+4;
                        basis.test.Tb = Tb;
                    case 'ipdg' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                        Pb = zeros(5,mesh.N);
                        Pb(1,:) = 1:2:2*mesh.N-1;  Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:)+1; Pb(4,:) = Pb(1,:)+2;
                        Pb(5,:) = Pb(1,:)+2;
                        basis.test.Pb = Pb0(reshape(Pb,1,5*mesh.N));
                        Tb = zeros(5,mesh.N);
                        Tb(1,:) = 1:5:5*mesh.N-4; Tb(2,:) = Tb(1,:)+2;
                        Tb(3,:) = Tb(1,:)+3;  Tb(4,:) = Tb(1,:)+1;
                        Tb(5,:) = Tb(1,:)+4;
                        basis.test.Tb = Tb;
                end
                 fclass = cell(7,5);
%                  fclass{1,1} = @(x) l1(x)^2*(l1(x)-l2(x))*(l1(x)+5*l2(x));
%                  fclass{1,2} = @(x) 16*l1(x)^2*l2(x)^2;
%                  fclass{1,3} = @(x) -l2(x)^2*(l1(x)-l2(x))*(5*l1(x)+l2(x));
%                  fclass{1,4} = @(x) 2*l1(x)^2*l2(x)*(l1(x)-l2(x));
%                  fclass{1,5} = @(x) 2*l1(x)*l2(x)^2*(l1(x)-l2(x));
                 fclass{1,1} = @(x) -1/4*(-1 + x).^2.*x.*(3 + 2*x);
                 fclass{1,2} = @(x) (-1 + x.^2).^2;
                 fclass{1,3} = @(x) -1/4*x.*(1 + x).^2.*(-3 + 2*x);
                 fclass{1,4} = @(x) -1/4*(-1 + x).^2.*x.*(1 + x);
                 fclass{1,5} = @(x) 1/4*(-1 + x).*x.*(1 + x).^2;
                 fclass{2,1} = @(x) 1/4*(-3 + 8*x).*(1 - x.^2);
                 fclass{2,2} = @(x) 4*x.*(-1 + x.^2);
                 fclass{2,3} = @(x) 1/4*(3 + 8*x).*(1 - x.^2);
                 fclass{2,4} = @(x) 1/4*(-1 + x).*(1 - x - 4*x.^2);
                 fclass{2,5} = @(x) 1/4*(1 + x).*(-1 - x + 4*x.^2);
                 fclass{3,1} = @(x) 2 + 3*x/2 - 6*x.^2;
                 fclass{3,2} = @(x) -4 + 12*x.^2;
                 fclass{3,3} = @(x) 2 - 3*x/2 - 6*x.^2;
                 fclass{3,4} = @(x) 1/2*(1 + 3*x - 6*x.^2);
                 fclass{3,5} = @(x) 1/2*(-1 + 3*x + 6*x.^2);
                 fclass{4,1} = @(x) 3/2 - 12*x;  
                 fclass{4,2} = @(x) 24*x;  
                 fclass{4,3} = @(x) -3/2*(1 + 8*x);
                 fclass{4,4} = @(x) 3/2 - 6*x;
                 fclass{4,5} = @(x) 3/2 + 6*x;
                 fclass{5,1} = @(x) -12;  
                 fclass{5,2} = @(x) 24;  
                 fclass{5,3} = @(x) -12;  
                 fclass{5,4} = @(x) -6; 
                 fclass{5,5} = @(x) 6;
                 fclass{6,1} = @(x) 0;  fclass{6,2} = @(x) 0;  
                 fclass{6,3} = @(x) 0;  fclass{6,4} = @(x) 0;
                 fclass{6,5} = @(x) 0;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;
                 basis.test.reffun = fclass;
            case 5  % Hermite fifth-order element: a1, da1, d^2a1, a2, da2, d^2a2
                switch femtype.type
                    case 'pgm' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                        Pb = zeros(3,mesh.N+1);
                        Pb(1,:) = 1:mesh.N+1; Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:);
                        Pb = reshape(Pb,1,3*mesh.N+3);
                        basis.test.Pb = Pb0(Pb);
                        Tb = zeros(6,mesh.N);
                        Tb(1,:) = 1:3:3*mesh.N-2; Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+1; Tb(4,:) = Tb(1,:)+4;
                        Tb(5,:) = Tb(1,:)+2; Tb(6,:) = Tb(1,:)+5;
                        basis.test.Tb = Tb;
                    case 'ipdg' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/mesh.N:mesh.right;
                        Pb = zeros(6,mesh.N);
                        Pb(1,:) = 1:mesh.N;  Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:); Pb(4,:) = Pb(1,:)+1;
                        Pb(5,:) = Pb(1,:)+1; Pb(6,:) = Pb(1,:)+1; 
                        basis.test.Pb = Pb0(reshape(Pb,1,6*mesh.N));
                        Tb = zeros(6,mesh.N);
                        Tb(1,:) = 1:6:6*mesh.N-5;  Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+1; Tb(4,:) = Tb(1,:)+4;
                        Tb(5,:) = Tb(1,:)+2; Tb(6,:) = Tb(1,:)+5;
                        basis.test.Tb = Tb;
                end
                 fclass = cell(7,6);
%                  fclass{1,1} = @(x) l1(x)^3*(l1(x)^2+5*l1(x)*l2(x)+10*l2(x)^2);
%                  fclass{1,2} = @(x) l2(x)^3*(10*l1(x)^2+5*l1(x)*l2(x)+l2(x)^2);
%                  fclass{1,3} = @(x) 2*l1(x)^3*l2(x)*(l1(x)+4*l2(x));
%                  fclass{1,4} = @(x) -2*l1(x)*l2(x)^3*(4*l1(x)+l2(x));
%                  fclass{1,5} = @(x) 2*l1(x)^3*l2(x)^2;
%                  fclass{1,6} = @(x) 2*l1(x)^2*l2(x)^3;
                 fclass{1,1} = @(x) -1/16*(-1 + x).^3.*(8 + 9*x + 3*x.^2);
                 fclass{1,2} = @(x) 1/16*(1 + x).^3.*(8 - 9*x + 3*x.^2);
                 fclass{1,3} = @(x) -1/16*(-1 + x).^3.*(1 + x).*(5 + 3*x);
                 fclass{1,4} = @(x) -1/16*(-1 + x).*(1 + x).^3.*(-5 + 3*x);
                 fclass{1,5} = @(x) -1/16*(-1 + x).^3.*(1 + x).^2;
                 fclass{1,6} = @(x) 1/16*(-1 + x).^2.*(1 + x).^3;
                 fclass{2,1} = @(x) -15/16*(-1 + x.^2).^2;
                 fclass{2,2} = @(x) 15/16*(-1 + x.^2).^2;
                 fclass{2,3} = @(x) -1/16*(-1 + x).^2.*(1 + 3*x).*(7 + 5*x);
                 fclass{2,4} = @(x) -1/16*(1 + x).^2.*(-1 + 3*x).*(-7 + 5*x);
                 fclass{2,5} = @(x) -1/16*(-1 + x).^2.*(1 + x).*(1 + 5*x);
                 fclass{2,6} = @(x) 1/16*(1 + x).^2.*(-1 + x).*(-1 + 5*x);
                 fclass{3,1} = @(x) -15/4*x.*(-1 + x.^2);
                 fclass{3,2} = @(x) 15/4*x.*(-1 + x.^2);
                 fclass{3,3} = @(x) -3/4*(-1 + 5*x).*(- 1 + x.^2);
                 fclass{3,4} = @(x) -3/4*(1 + 5*x).*(- 1 + x.^2);
                 fclass{3,5} = @(x) 1/4*(-1 + 3*x + 3*x.^2 - 5*x.^3);
                 fclass{3,6} = @(x) 1/4*(-1 - 3*x + 3*x.^2 + 5*x.^3);
                 fclass{4,1} = @(x) -15/4*(-1 + 3*x.^2); 
                 fclass{4,2} = @(x) 15/4*(-1 + 3*x.^2);  
                 fclass{4,3} = @(x) -3/4*(-5 - 2*x + 15*x.^2);
                 fclass{4,4} = @(x) -3/4*(-5 + 2*x + 15*x.^2);
                 fclass{4,5} = @(x) -3/4*(-1 - 2*x + 5*x.^2);
                 fclass{4,6} = @(x) 3/4*(-1 + 2*x + 5*x.^2); 
                 fclass{5,1} = @(x) -45*x/2;  
                 fclass{5,2} = @(x) 45*x/2;  
                 fclass{5,3} = @(x) 3/2 - 45*x/2;  
                 fclass{5,4} = @(x) -3/2*(1 + 15*x); 
                 fclass{5,5} = @(x) -3/2*(-1 + 5*x);
                 fclass{5,6} = @(x) 3/2*(1 + 5*x);
                 fclass{6,1} = @(x) -45/2;  
                 fclass{6,2} = @(x) 45/2;  
                 fclass{6,3} = @(x) -45/2;  
                 fclass{6,4} = @(x) -45/2;
                 fclass{6,5} = @(x) -15/2;  
                 fclass{6,6} = @(x) 15/2;
                 fclass{7,1} = @(x) 0;  fclass{7,2} = @(x) 0;  
                 fclass{7,3} = @(x) 0;  fclass{7,4} = @(x) 0;
                 fclass{7,5} = @(x) 0;  fclass{7,6} = @(x) 0;
                 basis.test.reffun = fclass;
            case 6  % Hermite sixth-order element: a1, da1, d^2a1, a2, da2, d^2a2, a3
                switch femtype.type
                    case 'pgm' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                        Pb = zeros(4,mesh.N+1);
                        Pb(1,:) = 1:2:2*mesh.N+1; Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:); Pb(4,:) = Pb(1,:)+1;
                        Pb = Pb(1:end-1);
                        basis.test.Pb = Pb0(reshape(Pb,1,4*mesh.N+3));
                        Tb = zeros(7,mesh.N);
                        Tb(1,:) = 1:4:4*mesh.N-3; Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+4; Tb(4,:) = Tb(1,:)+1;
                        Tb(5,:) = Tb(1,:)+5; Tb(6,:) = Tb(1,:)+2;
                        Tb(7,:) = Tb(1,:)+6;
                        basis.test.Tb = Tb;
                    case 'ipdg' 
                        Pb0 = mesh.left:(mesh.right-mesh.left)/(2*mesh.N):mesh.right;
                        Pb = zeros(7,mesh.N);
                        Pb(1,:) = 1:2:2*mesh.N-1;  Pb(2,:) = Pb(1,:);
                        Pb(3,:) = Pb(1,:); Pb(4,:) = Pb(1,:)+1;
                        Pb(5,:) = Pb(1,:)+2; Pb(6,:) = Pb(1,:)+2; 
                        Pb(7,:) = Pb(1,:)+2;
                        basis.test.Pb = Pb0(reshape(Pb,1,7*mesh.N));
                        Tb = zeros(7,mesh.N);
                        Tb(1,:) = 1:7:7*mesh.N-6;  Tb(2,:) = Tb(1,:)+3;
                        Tb(3,:) = Tb(1,:)+4; Tb(4,:) = Tb(1,:)+1;
                        Tb(5,:) = Tb(1,:)+5; Tb(6,:) = Tb(1,:)+2;
                        Tb(7,:) = Tb(1,:)+6;
                        basis.test.Tb = Tb;
                end
                 fclass = cell(7,7);
%                  fclass{1,1} = @(x) l1(x)^3*(l1(x)-l2(x))*(l1(x)^2+7*l1(x)*l2(x)+22*l2(x)^2);
%                  fclass{1,2} = @(x) 64*l1(x)^3*l2(x)^3;
%                  fclass{1,3} = @(x) -l2(x)^3*(l1(x)-l2(x))*(22*l1(x)^2+7*l1(x)*l2(x)+l2(x)^2);
%                  fclass{1,4} = @(x) l1(x)^3*l2(x)*(l1(x)-l2(x))*(2*l1(x)+12*l2(x));
%                  fclass{1,5} = @(x) l1(x)*l2(x)^3*(l1(x)-l2(x))*(12*l1(x)+2*l2(x));
%                  fclass{1,6} = @(x) 2*l1(x)^3*l2(x)^2*(l1(x)-l2(x));
%                  fclass{1,7} = @(x) -2*l1(x)^2*l2(x)^3*(l1(x)-l2(x));
                 fclass{1,1} = @(x) 1/16*(-1 + x).^3.*x.*(15 + 21*x + 8*x.^2);
                 fclass{1,2} = @(x) -(-1 + x.^2).^3;
                 fclass{1,3} = @(x) 1/16*x.*(1 + x).^3.*(15 - 21*x + 8*x.^2);
                 fclass{1,4} = @(x) 1/16*(-1 + x).^3.*x.*(1 + x).*(7 + 5*x);
                 fclass{1,5} = @(x) -1/16*(-1 + x).*x.*(1 + x).^3.*(-7 + 5*x);
                 fclass{1,6} = @(x) 1/16*(-1 + x).^3.*x.*(1 + x).^2;
                 fclass{1,7} = @(x) 1/16*(-1 + x).^2.*x.*(1 + x).^3;
                 fclass{2,1} = @(x) 3/16*(-5 + 16*x).*(-1 + x.^2).^2;
                 fclass{2,2} = @(x) -6*x.*(-1 + x.^2).^2;
                 fclass{2,3} = @(x) 3/16*(5 + 16*x).*(-1 + x.^2).^2;
                 fclass{2,4} = @(x) 1/16*(-1 + x).^2.*(-7 + 4*x + 45*x.^2 + 30*x.^3);
                 fclass{2,5} = @(x) -1/16*(1 + x).^2.*(7 + 4*x - 45*x.^2 + 30*x.^3);
                 fclass{2,6} = @(x) 1/16*(-1 + x).^2.*(-1 + 7*x.^2 + 6*x.^3);
                 fclass{2,7} = @(x) 1/16*(1 + x).^2.*(1 - 7*x.^2 + 6*x.^3);
                 fclass{3,1} = @(x) 3/4*(4 + 5*x - 24*x.^2 - 5*x.^3 + 20*x.^4);
                 fclass{3,2} = @(x) -6*(1 - 6*x.^2 + 5*x.^4);
                 fclass{3,3} = @(x) 3/4*(4 - 5*x - 24*x.^2 + 5*x.^3 + 20*x.^4);
                 fclass{3,4} = @(x) 3/8*(3 + 10*x - 28*x.^2 - 10*x.^3 + 25*x.^4);
                 fclass{3,5} = @(x) -3/8*(3 - 10*x - 28*x.^2 + 10*x.^3 + 25*x.^4);
                 fclass{3,6} = @(x) 1/8*(1 + 6*x - 12*x.^2 - 10*x.^3 + 15*x.^4);
                 fclass{3,7} = @(x) 1/8*(1 - 6*x - 12*x.^2 + 10*x.^3 + 15*x.^4);
                 fclass{4,1} = @(x) 3/4*(5 - 48*x - 15*x.^2 + 80*x.^3);  
                 fclass{4,2} = @(x) 24*x.*(3 - 5*x.^2);  
                 fclass{4,3} = @(x) 3/4*(-5 - 48*x + 15*x.^2 + 80*x.^3);
                 fclass{4,4} = @(x) 3/4*(5 - 28*x - 15*x.^2 + 50*x.^3);
                 fclass{4,5} = @(x) -3/4*(-5 - 28*x + 15*x.^2 + 50*x.^3);
                 fclass{4,6} = @(x) 3/4*(1 - 4*x - 5*x.^2 + 10*x.^3);
                 fclass{4,7} = @(x) 3/4*(-1 - 4*x + 5*x.^2 + 10*x.^3);
                 fclass{5,1} = @(x) 9/2*(-8 - 5*x + 40*x.^2);  
                 fclass{5,2} = @(x) 72 - 360*x.^2;  
                 fclass{5,3} = @(x) 9/2*(-8 + 5*x + 40*x.^2);  
                 fclass{5,4} = @(x) 3/2*(-14 - 15*x + 75*x.^2); 
                 fclass{5,5} = @(x) -3/2*(-14 + 15*x + 75*x.^2);
                 fclass{5,6} = @(x) 3/2*(-2 - 5*x + 15*x.^2);
                 fclass{5,7} = @(x) 3/2*(-2 + 5*x + 15*x.^2);
                 fclass{6,1} = @(x) -45/2 + 360*x;  
                 fclass{6,2} = @(x) -720*x;  
                 fclass{6,3} = @(x) 45/2 + 360*x;  
                 fclass{6,4} = @(x) -45/2 + 225*x;
                 fclass{6,5} = @(x) -45/2*(1 + 10*x);  
                 fclass{6,6} = @(x) -15/2 + 45*x;
                 fclass{6,7} = @(x) 15/2 + 45*x;
                 fclass{7,1} = @(x) 360;  
                 fclass{7,2} = @(x) -720;  
                 fclass{7,3} = @(x) 360;  
                 fclass{7,4} = @(x) 225;
                 fclass{7,5} = @(x) -450/2;  
                 fclass{7,6} = @(x) 45;
                 fclass{7,7} = @(x) 45;
                 basis.test.reffun = fclass;
        end
    basis.test.Nb = size(basis.test.Pb,2);
    basis.test.Nlb = size(basis.test.Tb,1);
end

end

newbasis = basis;

end
