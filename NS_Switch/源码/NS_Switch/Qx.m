function z=Qx(x,y,k,xk,yk,DET)

L1=((xk(2)-x)*(yk(3)-y)-(xk(3)-x)*(yk(2)-y))/DET;
L2=((x-xk(1))*(yk(3)-yk(1))-(xk(3)-xk(1))*(y-yk(1)))/DET;
L3=((xk(2)-xk(1))*(y-yk(1))-(x-xk(1))*(yk(2)-yk(1)))/DET;

if k==1
    z=(yk(2)-yk(3))*(4*L1-1)/DET;
elseif k==2
    z=(yk(3)-yk(1))*(4*L2-1)/DET;
elseif k==3
    z=(yk(1)-yk(2))*(4*L3-1)/DET;
elseif k==5
    z=4*((yk(3)-yk(1))*L3+(yk(1)-yk(2))*L2)/DET;
elseif k==6
    z=4*((yk(2)-yk(3))*L3+(yk(1)-yk(2))*L1)/DET;
elseif k==4
    z=4*((yk(2)-yk(3))*L2+(yk(3)-yk(1))*L1)/DET;
end
return
end