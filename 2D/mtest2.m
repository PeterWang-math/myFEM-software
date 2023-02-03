% -div(c·grad(u)) + a·u = f ;
% \omega = (-1,1)x(-1,1) 
% Dirichlet boundary condition
% exact solution: u(x,y)=exp(x+y)
g = [ 2     2     2     2
     -1     1     1    -1
      1     1    -1    -1
      1     1    -1    -1
      1    -1    -1     1
      0     0     0     0
      1     1     1     1 ];
b = [1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     8     8     8     9
    48    48    48    48
    48    48    48    48
    49    49    49    49
   101   101   101   101
   120   120   120   120
   112   112   112   112
    40    40    40    40
   120    49   120    45
    43    43    45    49
    49   121    49    43
    41    41    41   121
     0     0     0    41 ];
c='1'; a='0'; f= '-2*exp(x+y)';
n=10;
[p,e,t]=poimesh(g,n);
pdemesh(p,e,t);
[K,~,F,~,~,~,~]=assempde(b,p,e,t,c,a,f);
u=assempde(b,p,e,t,c,a,f);
pdesurf(p,t,u);
