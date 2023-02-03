function u = lshape(x,y)
r = sqrt(x.^2+y.^2);
theta = atan2(y,x);
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
u = r.^(2/3).*sin(2*theta/3);
end

