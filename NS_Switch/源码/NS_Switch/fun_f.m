function z=fun_f(nu,t,x,y)
n=length(x);
%z1=(-2*nu*x.^2-2*nu*y.^2-nu*exp(-y)+pi^2*cos(pi*x).*cos(2*pi*y)).*cos(2*pi*t)-2*pi*(x.^2.*y.^2+exp(-y)).*sin(2*pi*t);  
%z1=z1+((x.^2.*y.^2+exp(-y))*2.*x.*y.^2+(-2*x.*y.^3/3+2-pi*sin(pi*x)).*(2*x.^2.*y-exp(-y))).*cos(2*pi*t).*cos(2*pi*t);
%z2=(4*nu*x.*y-nu*pi^3*sin(pi*x)+2*pi.*(2-pi*sin(pi*x)).*sin(2*pi*y)).*cos(2*pi*t)-2*pi*(-2*x.*y.^3/3+2-pi*sin(pi*x)).*sin(2*pi*t);   
%z2=z2+((x.^2.*y.^2+exp(-y)).*(-2*y.^3/3-pi^2*cos(pi*x))+(-2*x.*y.^3/3+2-pi*sin(pi*x)).*(-2*x.*y.^2)).*cos(2*pi*t).*cos(2*pi*t);
%z1=2*t*sin(3*pi*x).*cos(2*pi*y);
%z2=2*t*cos(3*pi*x).*sin(2*pi*y);
z1=zeros(n,1);
z2=zeros(n,1);
z=[z1;z2];
return
end