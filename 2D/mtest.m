% -div(c·grad(u)) + a·u = f ;
% \omega = (-1,1)x(-1,1) 
% Dirichlet zero boundary condition
g = [2     2     2     2
    -1     1     1    -1
     1     1    -1    -1
    -1    -1     1     1
    -1     1     1    -1
     1     1     1     1
     0     0     0     0];
b = [1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
    48    48    48    48
    48    48    48    48
    49    49    49    49
    48    48    48    48 ];
c='(2+1.8*sin(2*pi*x/1e-1))./(2+1.8*cos(2*pi*y/1e-1))+(2+1.8*sin(2*pi*y/1e-1))./(2+1.8*sin(2*pi*x/1e-1))';  
%c='(2+1.8*sin(2*pi*x/1e-2))./(2+1.8*cos(2*pi*y/1e-2))+(2+1.8*sin(2*pi*y/1e-2))./(2+1.8*sin(2*pi*x/1e-2))';
a='0';
f='1';
n=1000;
[p,e,t]=poimesh(g,n);
% figure(1);
% pdemesh(p,e,t);
[K,~,F,~,~,~,~]=assempde(b,p,e,t,c,a,f);
u=assempde(b,p,e,t,c,a,f);
figure(2);
pdesurf(p,t,u);

% L-shape
% u=sqrt(sum(x.^2+y.^2,2)).^(2/3).*sin(2*(atan2(y,x)>=0).*atan2(y,x) + (atan2(y,x)<0).*(atan2(y,x)+2*pi)/3);
% %x=0
% sqrt(sum(y.^2,2)).^(2/3).*sin(2*(atan2(y,0)>=0).*atan2(y,0) + (atan2(y,0)<0).*(atan2(y,0)+2*pi)/3);
% %x=1
% sqrt(sum(1+y.^2,2)).^(2/3).*sin(2*(atan2(y,1)>=0).*atan2(y,1) + (atan2(y,1)<0).*(atan2(y,1)+2*pi)/3);
% %x=-1
% sqrt(sum(1.^2+y.^2,2)).^(2/3).*sin(2*(atan2(y,-1)>=0).*atan2(y,-1) + (atan2(y,-1)<0).*(atan2(y,-1)+2*pi)/3);
% %y=0
% sqrt(sum(x.^2,2)).^(2/3).*sin(2*(atan2(0,x)>=0).*atan2(0,x) + (atan2(0,x)<0).*(atan2(0,x)+2*pi)/3);
% %y=1
% sqrt(sum(x.^2+1,2)).^(2/3).*sin(2*(atan2(1,x)>=0).*atan2(1,x) + (atan2(1,x)<0).*(atan2(1,x)+2*pi)/3);
% %y=-1
% sqrt(sum(x.^2+1,2)).^(2/3).*sin(2*(atan2(-1,x)>=0).*atan2(-1,x) + (atan2(-1,x)<0).*(atan2(-1,x)+2*pi)/3);



