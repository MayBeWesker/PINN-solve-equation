function z=u2_bry(t,x,y,boundarynodes)

% boundaryU2.gam_1=@(t,x,y) (-2*x.*y.^3/3+2-pi*sin(pi*x)).*cos(2*pi*t);
% boundaryU2.gam_2=@(t,x,y) (-2*x.*y.^3/3+2-pi*sin(pi*x)).*cos(2*pi*t);
% boundaryU2.gam_3=@(t,x,y) (-2*x.*y.^3/3+2-pi*sin(pi*x)).*cos(2*pi*t);
% boundaryU2.gam_4=@(t,x,y) (-2*x.*y.^3/3+2-pi*sin(pi*x)).*cos(2*pi*t);
boundaryU2.gam_1=@(t,x,y) 0;
boundaryU2.gam_2=@(t,x,y) 0;
boundaryU2.gam_3=@(t,x,y) 0;
boundaryU2.gam_4=@(t,x,y) 0;

z=zeros(length(x),1);

gam_1_index=find(boundarynodes(2,:)==1);
z(gam_1_index)=boundaryU2.gam_1(t,x(gam_1_index),y(gam_1_index));

gam_2_index=find(boundarynodes(2,:)==2);
z(gam_2_index)=boundaryU2.gam_2(t,x(gam_2_index),y(gam_2_index));

gam_3_index=find(boundarynodes(2,:)==3);
z(gam_3_index)=boundaryU2.gam_3(t,x(gam_3_index),y(gam_3_index));

gam_4_index=find(boundarynodes(2,:)==4);
z(gam_4_index)=boundaryU2.gam_4(t,x(gam_4_index),y(gam_4_index));

return
end