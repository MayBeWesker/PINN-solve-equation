function z=boundary_pts(P,x,y)
% z(1,k):=boundary_pts(1,k) is the type of the kth boundary point
% Dirichlet (-1), Neumann(-2), Robin(-3)
% z(2,k):=boundary_pts(2,k) is the index of which boundary does the kth point belong to
% z=0: internal nodes
z=zeros(2,size(P,2));
    
gam_1=find(P(1,:)==x(1));
z(1,gam_1)=-1;  % Dirichlet
z(2,gam_1)=1;  % 第一条边
    
gam_2=find(P(2,:)==y(1));
z(1,gam_2)=-1;  % Dirichlet
z(2,gam_2)=2;  % 第二条边
    
gam_3=find(P(1,:)==x(2));
z(1,gam_3)=-1;  % Dirichlet
z(2,gam_3)=3;  % 第三条边
    
gam_4=find(P(2,:)==y(2));
z(1,gam_4)=-1;  % Dirichlet
z(2,gam_4)=4;  % 第四条边
     
return
end