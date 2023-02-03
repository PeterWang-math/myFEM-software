function u = lshapeddyy(x,y)
r = sqrt(x.^2+y.^2);
theta = atan2(y,x);
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
u = r.^(-10/3)/9.*( -4*x.*y.*cos(2*theta/3) + 2*(x.^2-y.^2).*sin(2*theta/3) );
end

