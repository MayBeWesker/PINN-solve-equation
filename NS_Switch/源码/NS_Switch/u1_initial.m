function z=u1_initial(x,y)
%z=x.^2.*y.^2+exp(-y);
z=(-10*x.*x.*(1-x).*(1-x)).*(20*y.*(1-y).*(1-y)-20*y.*y.*(1-y));
return
end
