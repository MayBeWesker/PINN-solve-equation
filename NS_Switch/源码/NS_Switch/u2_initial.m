function z=u2_initial(x,y)
%z=(-2*x.*y.^3/3+2-pi*sin(pi*x));
z=(20*x.*(1-x).*(1-x)-20*x.*x.*(1-x)).*(10*y.*y.*(1-y).*(1-y));
return
end