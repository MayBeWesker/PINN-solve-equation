function [u1,u2,P_n]=ns_pde_solver(theta,T,m,nu,u1_initial,u2_initial,p_initial,...
                                    Pb,Tb,Nb,N1,elm_bry)
num.nodes=length(Pb(1,:)); % the number of finite element nodes
num.elem=length(Tb(1,:));  % the number of finite elements
dt=T/m; % time step
width=max(Pb(1,:));   % x_max
Height=max(Pb(2,:));  % y_max

domain.x=[0,width]; % type of nodes
domain.y=[0,Height];
% find bry points about special boundary and type of the boundary pts
boundarynodes=boundary_pts(Pb,domain.x,domain.y);

%---------------------------------
%linear term form
W1=[0.225,0.1323941527*ones(1,3),0.1259391805*ones(1,3)]; 
alpha=[0.0597158717, 0.7974269853];
beta =[0.4701420641, 0.1012865073];

P1=[1/3,1/3,1/3;
    alpha(1),beta(1),beta(1);
    beta(1),alpha(1),beta(1);
    beta(1),beta(1),alpha(1);
    alpha(2),beta(2),beta(2);
    beta(2),alpha(2),beta(2);
    beta(2),beta(2),alpha(2)];

DT = zeros(num.elem,1);   % 纪录各单元面积
Index.nodes=zeros(1,6*num.elem); Index.bry=zeros(1,num.elem);
Gbl.Pijkx=spalloc(6*num.elem,36*num.elem,6*36*num.elem);
Gbl.Pijky=spalloc(6*num.elem,36*num.elem,6*36*num.elem);
Gbl.Tijkx=spalloc(6*num.elem,36*num.elem,6*36*num.elem);
Gbl.Tijky=spalloc(6*num.elem,36*num.elem,6*36*num.elem);

Val_Dx_Q=zeros(6,3,num.elem);   % 纪录6个基对x一阶导在3个高斯节点（各边中点）的值
Val_Dy_Q=zeros(6,3,num.elem);   % 纪录6个基对y一阶导在3个高斯节点（各边中点）的值
Val_Q=zeros(6,7,num.elem);      % 纪录6个基在7个高斯节点的值
Gbl_Dx_Q=zeros(6,7,num.elem);   % 纪录6个基对x一阶导在7个高斯节点的值
Gbl_Dy_Q=zeros(6,7,num.elem);   % 纪录6个基对y一阶导在7个高斯节点的值
Val_L=zeros(3,3,num.elem);      % 记录线性基函数在3个高斯节点的值

for i=1:num.elem
    xk=Pb(1,Tb(:,i)); yk=Pb(2,Tb(:,i));
    xx=P1*xk(1:3)';   yy=P1*yk(1:3)';
    AA=[ones(3,1),xk(1:3)',yk(1:3)'];
    DET=det(AA);  % Area of the triangle
    DT(i)=DET;
    
    Index.nodes((i-1)*6+1:i*6)=Tb(:,i);
     % elevaluate the TH basis(2-order basis) and derivatives at the nodes firstly
    for j=1:6
        Val_Dx_Q(j,:,i) = [ Qx(xk(4),yk(4),j,xk(1:3),yk(1:3),DET),...  
                            Qx(xk(5),yk(5),j,xk(1:3),yk(1:3),DET),...
                            Qx(xk(6),yk(6),j,xk(1:3),yk(1:3),DET)];  % D_x(Q_j) at nodes__{i,4},nodes__{i,5},nodes_{i,6}
           
        Val_Dy_Q(j,:,i) = [ Qy(xk(4),yk(4),j,xk(1:3),yk(1:3),DET),...
                            Qy(xk(5),yk(5),j,xk(1:3),yk(1:3),DET),...
                            Qy(xk(6),yk(6),j,xk(1:3),yk(1:3),DET)];  % D_y(Q_j) at nodes__{i,4},nodes__{i,5},nodes_{i,6}
        
        Val_Q(j,:,i)    = [ Q(xx(1),yy(1),j,xk(1:3),yk(1:3),DET),...
                            Q(xx(2),yy(2),j,xk(1:3),yk(1:3),DET),...
                            Q(xx(3),yy(3),j,xk(1:3),yk(1:3),DET),...
                            Q(xx(4),yy(4),j,xk(1:3),yk(1:3),DET),...
                            Q(xx(5),yy(5),j,xk(1:3),yk(1:3),DET),...
                            Q(xx(6),yy(6),j,xk(1:3),yk(1:3),DET),...
                            Q(xx(7),yy(7),j,xk(1:3),yk(1:3),DET)];   % Q_j at nodes (xx_i,yy_i) 
                        
         Gbl_Dx_Q(j,:,i)= [ Qx(xx(1),yy(1),j,xk(1:3),yk(1:3),DET),...
                            Qx(xx(2),yy(2),j,xk(1:3),yk(1:3),DET),...
                            Qx(xx(3),yy(3),j,xk(1:3),yk(1:3),DET),...
                            Qx(xx(4),yy(4),j,xk(1:3),yk(1:3),DET),...
                            Qx(xx(5),yy(5),j,xk(1:3),yk(1:3),DET),...
                            Qx(xx(6),yy(6),j,xk(1:3),yk(1:3),DET),...
                            Qx(xx(7),yy(7),j,xk(1:3),yk(1:3),DET)];  % DxQ_j at nodes (xx_i,yy_i)
                        
         Gbl_Dy_Q(j,:,i)= [ Qy(xx(1),yy(1),j,xk(1:3),yk(1:3),DET),...
                            Qy(xx(2),yy(2),j,xk(1:3),yk(1:3),DET),...
                            Qy(xx(3),yy(3),j,xk(1:3),yk(1:3),DET),...
                            Qy(xx(4),yy(4),j,xk(1:3),yk(1:3),DET),...
                            Qy(xx(5),yy(5),j,xk(1:3),yk(1:3),DET),...
                            Qy(xx(6),yy(6),j,xk(1:3),yk(1:3),DET),...
                            Qy(xx(7),yy(7),j,xk(1:3),yk(1:3),DET)];  % DyQ_j at nodes (xx_i,yy_i)     
         if j<=3 
            Val_L(j,:,i)= [ L(xk(4),yk(4),j,xk(1:3),yk(1:3),DET),... % linear basis 
                            L(xk(5),yk(5),j,xk(1:3),yk(1:3),DET),...
                            L(xk(6),yk(6),j,xk(1:3),yk(1:3),DET)];
         end
    end
    for l1=1:6
       for l2=1:6
        % (phi_j,(phi_k)_x, phi_i)
        Gbl.Pijkx((i-1)*6+1:i*6,(i-1)*36+l2+(l1-1)*6)=((repmat(Val_Q(l1,:,i).*Val_Q(l2,:,i),6,1).*Gbl_Dx_Q(:,:,i)))*W1'*DT(i)/2;
        % (phi_j,(phi_k)_y, phi_i) 
        Gbl.Pijky((i-1)*6+1:i*6,(i-1)*36+l2+(l1-1)*6)=((repmat(Val_Q(l1,:,i).*Val_Q(l2,:,i),6,1).*Gbl_Dy_Q(:,:,i)))*W1'*DT(i)/2;
        % (phi_k,(phi_j)_x, phi_i)
        Gbl.Tijkx((i-1)*6+1:i*6,(i-1)*36+l2+(l1-1)*6)=((repmat(Val_Q(l1,:,i).*Gbl_Dx_Q(l2,:,i),6,1).*Val_Q(:,:,i)))*W1'*DT(i)/2;
        % (phi_k,(phi_j)_y, phi_i)
        Gbl.Tijky((i-1)*6+1:i*6,(i-1)*36+l2+(l1-1)*6)=((repmat(Val_Q(l1,:,i).*Gbl_Dy_Q(l2,:,i),6,1).*Val_Q(:,:,i)))*W1'*DT(i)/2; 
       end
   end     
end

% assemble matrix
K1=spalloc(num.nodes,num.nodes,num.nodes); Me=K1;
H1=zeros(6,6,num.elem); H2=H1;
K3=zeros(num.nodes,N1); K4=K3;
for i=1:num.elem
    for j=1:6
        % K1 = a(phi_j,phi_i)
        H1(j,:,i)=nu*(sum(repmat(Val_Dx_Q(j,:,i),6,1).*Val_Dx_Q(:,:,i),2)...
                     +sum(repmat(Val_Dy_Q(j,:,i),6,1).*Val_Dy_Q(:,:,i),2))'*DT(i)/6; 
        % M = (phi_j,phi_i)
        H2(j,:,i) = W1*(repmat(Val_Q(j,:,i),6,1).*Val_Q(:,:,i))'*DT(i)/2;  
        % K3 =b（varphi_i,(phi_j)_x）= K5'
        KK3(j,:)=sum(repmat(Val_Dx_Q(j,:,i),3,1).*Val_L(:,:,i),2)'*DT(i)/6; 
        % K4 =b（varphi_i,(phi_j)_y）= K6'
        KK4(j,:)=sum(repmat(Val_Dy_Q(j,:,i),3,1).*Val_L(:,:,i),2)'*DT(i)/6;
    end
    K1(Tb(:,i),Tb(:,i))=K1(Tb(:,i),Tb(:,i))+H1(:,:,i);
    Me(Tb(:,i),Tb(:,i))=Me(Tb(:,i),Tb(:,i))+H2(:,:,i);
    K3(Tb(:,i),Tb(1:3,i))=K3(Tb(:,i),Tb(1:3,i))+KK3;
    K4(Tb(:,i),Tb(1:3,i))=K4(Tb(:,i),Tb(1:3,i))+KK4;
end
% the linear terms are written as a (2*num.nodes+N1)*(2*num.nodes+N1) matrix
%     [K1   0  -K3]  [M/(theta*dt)       0       0]  [K1   0   -K3]  [M/(theta*dt)       0       0]  
% K = [0    K1 -K4]+ [    0        M/(theta*dt)  0]= [0    K1  -K4]+ [    0        M/(theta*dt)  0]
%     [-K5 -K6  0 ]  [    0              0       0]  [-K3' -K4'  0]  [    0              0       0]
O1=zeros(num.nodes,num.nodes);
O2=zeros(num.nodes,N1);
O3=zeros(N1,N1);
KK=[K1 O1 -K3; O1 K1 -K4; K3' K4' O3];
M=[Me O1 O2 ; O1 Me  O2;  O2' O2' O3];
K=M/theta/dt+KK;

%-------------------------------------------------------
% deal with the Dirichlet boundary condition
ind_bry = find(boundarynodes(1,:)==-1);  % Dirichlet边界点的global index
%----------------------------------------

%----------------------------------------
% initial conditions t=0
u1=u1_initial;
u2=u2_initial;
P_n=p_initial;
%P_n(1:N1,1)=p_initial(Pb(1,1:N1)',Pb(2,1:N1)');
% Boundary_U1=u1_bry(0,Pb(1,:),Pb(2,:),boundarynodes);
% u1(ind_bry,1)=Boundary_U1(ind_bry);
% Boundary_U2=u2_bry(0,Pb(1,:),Pb(2,:),boundarynodes);
% u2(ind_bry,1)=Boundary_U2(ind_bry);

%----------------------------------------
%nonlinear term and solve the ns eqaution
for ti=1:m
    %ti
    Gbl.un=u1(Index.nodes,ti);  
    Gbl.vn=u2(Index.nodes,ti);
    
    Gbl.K11=Gbl.un'*Gbl.Pijkx; % 1*(36*num.elem)
    Gbl.C11=reshape(Gbl.K11',6,6*num.elem); % 6*(6*num.elem) % KN1
    
    Gbl.K12=Gbl.un'*Gbl.Pijky; % 1*(36*num.elem) 
    Gbl.C12=reshape(Gbl.K12',6,6*num.elem); % 6*(6*num.elem) % KN3
    
    Gbl.K21=Gbl.vn'*Gbl.Pijkx; % 1*(36*num.elem)
    Gbl.C21=reshape(Gbl.K21',6,6*num.elem); % 6*(6*num.elem) % KN2
    
    Gbl.K22=Gbl.vn'*Gbl.Pijky; % 1*(36*num.elem)
    Gbl.C22=reshape(Gbl.K22',6,6*num.elem); % 6*(6*num.elem) % KN4
    
    Gbl.KC1=Gbl.un'*Gbl.Tijkx; %1*(36*num.elem)
    Gbl.CC1=reshape(Gbl.KC1',6,6*num.elem);  % 6*(6*num.elem) % KN5
    
    Gbl.KC2=Gbl.vn'*Gbl.Tijky; %1*(36*num.elem)
    Gbl.CC2=reshape(Gbl.KC2',6,6*num.elem);  % 6*(6*num.elem) % KN6
    
    KN1=zeros(num.nodes,num.nodes);
    KN2=zeros(num.nodes,num.nodes);
    KN3=zeros(num.nodes,num.nodes);
    KN4=zeros(num.nodes,num.nodes);
    KN5=zeros(num.nodes,num.nodes);
    KN6=zeros(num.nodes,num.nodes);
    
    % curl=(grad \times u，grad \times u)
    CURL=0; 
    
    % Assembe matrix for nonlinear term
    for i=1:num.elem
        C11=Gbl.C11(:,(i-1)*6+1:i*6)';   % 6*6  % KN1的每个单元上的6*6矩阵
        C12=Gbl.C12(:,(i-1)*6+1:i*6)';   % 6*6  % KN3的每个单元上的6*6矩阵
        C21=Gbl.C21(:,(i-1)*6+1:i*6)';   % 6*6  % KN2的每个单元上的6*6矩阵
        C22=Gbl.C22(:,(i-1)*6+1:i*6)';   % 6*6  % KN4的每个单元上的6*6矩阵
        CC1=(Gbl.CC1(:,(i-1)*6+1:i*6))'; % 6*6  % KN5的每个单元上的6*6矩阵
        CC2=(Gbl.CC2(:,(i-1)*6+1:i*6))'; % 6*6  % KN6的每个单元上的6*6矩阵
        
        KN1(Tb(:,i),Tb(:,i))=KN1(Tb(:,i),Tb(:,i))+C11;
        KN2(Tb(:,i),Tb(:,i))=KN2(Tb(:,i),Tb(:,i))+C21;
        KN3(Tb(:,i),Tb(:,i))=KN3(Tb(:,i),Tb(:,i))+C12;
        KN4(Tb(:,i),Tb(:,i))=KN4(Tb(:,i),Tb(:,i))+C22;
        KN5(Tb(:,i),Tb(:,i))=KN5(Tb(:,i),Tb(:,i))+CC1;
        KN6(Tb(:,i),Tb(:,i))=KN6(Tb(:,i),Tb(:,i))+CC2;
        
        CURL=CURL+sum((u2(Tb(:,i),ti)'*Val_Dx_Q(:,:,i)-u1(Tb(:,i),ti)'*Val_Dy_Q(:,:,i)).^2)*DT(i)/6; 
    end
    
    % assemble matrice
    % the nonlinear terms are written as a (2*num.nodes+N1)*(2*num.nodes+N1) matrix
    %     [KN1+KN5+KN6      KN3      0]
    % KN =[    KN2      KN4+KN5+KN6  0]
    %     [     0            0       0]
    
    KN=[KN1+KN5+KN6 KN3 O2; KN2 KN4+KN5+KN6 O2; O2' O2' O3];
%     KN(ind_bry,:)=0;
%     KN(ind_bry+num.nodes,:)=0;
    A=K+KN;
    A(ind_bry,:)=0; 
    A(ind_bry+num.nodes,:)=0;
    for j=1:length(ind_bry)
        A(ind_bry(j),ind_bry(j))=1;
        A(num.nodes+ind_bry(j),num.nodes+ind_bry(j))=1;
    end
    
    %right-hand side
    %b=F+BNX^(m-1)+M/(theta*dt)X^(m-1)
    %     [KN5+KN6  0     0][U1]
    % BN =[0     KN5+KN6  0][U2]
    %     [0        0     0][P_n]
    bN1=zeros(num.nodes,1); 
    bN2=zeros(num.nodes,1);
    bN3=zeros(num.nodes,1);
    bN4=zeros(num.nodes,1);
    
    bN1=KN5*u1(:,ti);
    bN2=KN6*u1(:,ti);
    bN3=KN5*u2(:,ti);
    bN4=KN6*u2(:,ti);
    BN=[bN1+bN2; bN3+bN4; zeros(N1,1)];
    
    % force term
    F_theta=theta*fun_f(nu,ti*dt,Pb(1,:)',Pb(2,:)')...
           +(1-theta)*fun_f(nu,(ti-1)*dt,Pb(1,:)',Pb(2,:)');
    F=[F_theta;zeros(N1,1)];
    R=M/theta/dt*[u1(:,ti); u2(:,ti); P_n(:,ti)]+M*F+BN;
    % 右端项边界处理
    Boundary_U1=(1-theta).*u1_bry((ti-1)*dt,Pb(1,:),Pb(2,:),boundarynodes)...
               +theta*u1_bry(ti*dt,Pb(1,:),Pb(2,:),boundarynodes);
    R(ind_bry,1)=Boundary_U1(ind_bry);%对u1所对的R边界处理
    Boundary_U2=(1-theta)*u2_bry((ti-1)*dt,Pb(1,:),Pb(2,:),boundarynodes)...
               +theta*u2_bry(ti*dt,Pb(1,:),Pb(2,:),boundarynodes);
    R(num.nodes+ind_bry,1)=Boundary_U2(ind_bry);%对u2所对的R边界处理
    
    
    % deal with b(u^(m+1)_(theta),q)=0 if m>0; 
    %           b(u^(m+1)_(theta),q)=(1-theta)b(u^(0),q)) if m=0;
    if ti==1
        R(2*num.nodes+1:end)=(1-theta).*[K3' K4']*[u1(:,ti); u2(:,ti)];
    end
    
     % evaluate p at one given pt
    A(2*num.nodes+ind_bry(1),:)=0;
    A(2*num.nodes+ind_bry(1),2*num.nodes+ind_bry(1))=1; 
    R(2*num.nodes+ind_bry(1))=0;
    %theta*p_exact_solution(ti*dt,Pb(1,ind_bry(1)),Pb(2,ind_bry(1)))...
    %                         +(1-theta)*p_exact_solution((ti-1)*dt,Pb(1,ind_bry(1)),Pb(2,ind_bry(1)));
  
    w=sparse(A)\R;  % solve for (u,v,p)
    u1(:,ti+1) = (w(1:num.nodes) - (1-theta)*u1(:,ti))/theta;
    u2(:,ti+1) = (w(num.nodes+1:2*num.nodes)- (1-theta)*u2(:,ti))/theta;
    P_n(:,ti+1)= (w(2*num.nodes+1:end)-(1-theta)*P_n(:,ti))/theta;
   
end


return
end