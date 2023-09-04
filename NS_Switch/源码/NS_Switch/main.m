function main %保存的实验中时间步长1/128，空间步长1/64
nu = [0.01;0.02];     % viscosity雷诺数在全局时间上的变化
t=[0,0.5;0.5,1];   % 单位时间简单分成两段
m0 = 16; m=m0/length(nu);  % time intervals 时间步长
left=0; right=1; bottom=0; top=1; % space intervals 空间domain
h_partition=[1/16,1/16]; %空间步长
N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
nbn=((N1_partition+1)*(N2_partition+1));
theta=1/2;

[PP,NN,Nb,N1,elm_bry,e]=domain(left,right,bottom,top,h_partition);
Pb=PP;  Tb=NN; num.nodes = length(PP);
u1(1:num.nodes,1)=u1_initial(Pb(1,:)',Pb(2,:)');
u2(1:num.nodes,1)=u2_initial(Pb(1,:)',Pb(2,:)');
p(1:N1,1)=p_initial(Pb(1,1:N1)',Pb(2,1:N1)');

U0 = zeros(length(Nb),m0+1); V0 = zeros(length(Nb),m0+1);
U = zeros(N1_partition+1,N1_partition+1,2*(m+1)-1);
V = zeros(N1_partition+1,N1_partition+1,2*(m+1)-1);
P = zeros(N1_partition+1,N1_partition+1,2*(m+1)-1);
for i = 1:length(nu)
    [u1,u2,p] = Navier_stokes(theta,t(:,i),u1(:,end),u2(:,end),p(:,end),m,nu(i),Pb,Tb,Nb,N1,elm_bry);
    step=m+1;
    if i ~= 1
            u1 = u1(:,2:end);u2 = u2(:,2:end);p = p(:,2:end);
            step=m;U0(:,m*(i-1)+2:m*i+1) = u1;V0(:,m*(i-1)+2:m*i+1) = u2;
    else
            U0(:,1:step)=u1;V0(:,1:step)=u2; P0(:,1:step)=p;
    end

    for j = 1: step
        U(:,:,(i-1)*(m+1)+j) = reshape(u1(1:nbn,j),[N1_partition+1 N2_partition+1]);
        V(:,:,(i-1)*(m+1)+j) = reshape(u2(1:nbn,j),[N1_partition+1 N2_partition+1]);
        P(:,:,(i-1)*(m+1)+j) = reshape(p(1:nbn,j),[N1_partition+1 N2_partition+1]);
    end
    
    %绘制0.5s 和 1s的图
    figure
    quiver(Pb(1,1:N1)',Pb(2,1:N1)',u1(1:N1,end),u2(1:N1,end),'color','b','linewidth',1.5);
    axis tight
    axis equal
    
    figure
    pdeplot(Pb,e,Tb,'xydata',u1(:,end).^2+u2(:,end).^2,'zdata',u1(:,end).^2+u2(:,end).^2,'colormap','jet','mesh','off');view(2)
    axis tight
    axis equal
end
%保存在格子点上的数据U V P , U0 V0 P0是有限元划分点原始数据
save udata U V P U0 V0 P0 
end
