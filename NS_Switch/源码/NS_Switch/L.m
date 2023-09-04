function z=L(x,y,j,xk,yk,DET)

if j==1
    z=((xk(2)-x)*(yk(3)-y)-(xk(3)-x)*(yk(2)-y))/DET;
elseif j==2
    z=((x-xk(1))*(yk(3)-yk(1))-(xk(3)-xk(1))*(y-yk(1)))/DET;
elseif j==3
    z=((xk(2)-xk(1))*(y-yk(1))-(x-xk(1))*(yk(2)-yk(1)))/DET;
end
return
end