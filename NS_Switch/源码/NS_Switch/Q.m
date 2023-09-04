function z=Q(x,y,k,xk,yk,DET)

L1=((xk(2)-x)*(yk(3)-y)-(xk(3)-x)*(yk(2)-y))/DET;
L2=((x-xk(1))*(yk(3)-yk(1))-(xk(3)-xk(1))*(y-yk(1)))/DET;
L3=((xk(2)-xk(1))*(y-yk(1))-(x-xk(1))*(yk(2)-yk(1)))/DET;

if k==1 z=L1*(2*L1-1);
elseif k==2  z=L2*(2*L2-1);
elseif k==3  z=L3*(2*L3-1);
elseif k==4  z=4*L1*L2;
elseif k==5  z=4*L2*L3;
elseif k==6  z=4*L1*L3;
end
return
end