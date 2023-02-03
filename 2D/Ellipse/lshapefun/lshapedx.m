function u = lshapedx(x,y)
r = sqrt(x.^2+y.^2);
theta = atan2(y,x);
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
u = 2/3*r.^(-1/3).*sin(2*theta/3).*x./r ...
    - 2/3*r.^(2/3).*cos(2*theta/3).*y./r.^2;
end

