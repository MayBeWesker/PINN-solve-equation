function [PP,NN,Nb,N1,elm_bry,e]=domain(left,right,bottom,top,h_partition)
%load SquareDomain p e t
%load Lshape
[p,e,t]=generatePET_triangle(left,right,bottom,top,h_partition);
p=p; % the finite element nodes 的坐标
e=e(3:4,:); % boundary edge index
t=t; %每个单元上所有节点的编号阵

num.elem=length(t(1,:));
Nb=zeros(1,length(p));
Nb(e(1,:))=1; Nb(e(2,:))=1; 
% Nb：边界点的标识阵 边界点Nb(i)记为1，否则记为0

N1=length(p);

  i=1;
  s4=(p(:,t(1,i))+p(:,t(2,i)))/2;
  s5=(p(:,t(2,i))+p(:,t(3,i)))/2;
  s6=(p(:,t(3,i))+p(:,t(1,i)))/2;
  s=[s4,s5,s6];
  P=[N1+1,N1+2,N1+3]';
  N=N1+3;
  
  if Nb(t(1,i))==1 && Nb(t(2,i))==1
      Nb(N1+1)=1;
  else
      Nb(N1+1)=0;
  end
  
  if Nb(t(2,i))==1 && Nb(t(3,i))==1
      Nb(N1+2)=1;
  else
      Nb(N1+2)=0;
  end
  
  if Nb(t(3,i))==1 && Nb(t(1,i))==1
      Nb(N1+3)=1;
  else
      Nb(N1+3)=0;
  end

  for i=2:num.elem
      s4=(p(:,t(1,i))+p(:,t(2,i)))/2;
      s5=(p(:,t(2,i))+p(:,t(3,i)))/2;
      s6=(p(:,t(3,i))+p(:,t(1,i)))/2;
  
      I4=find_col(s,s4);
      if I4==0
          s=[s,s4]; P1=N+1; N=N+1;
          if Nb(t(1,i))==1 && Nb(t(2,i))==1
              Nb(N)=1;
          else
              Nb(N)=0;
          end
      else
          s=s; P1=N1+I4(2);
      end
      
      I5=find_col(s,s5);
      if I5==0
          s=[s,s5]; P2=N+1;N=N+1;
          if Nb(t(2,i))==1 && Nb(t(3,i))==1
              Nb(N)=1;
          else
              Nb(N)=0;
          end
      else
          s=s; P2=N1+I5(2);
      end
      
      I6=find_col(s,s6);
      if I6==0
          s=[s,s6]; P3=N+1;N=N+1;
          if Nb(t(3,i))==1 && Nb(t(1,i))==1
              Nb(N)=1;
          else
              Nb(N)=0;
          end
      else
          s=s; P3=N1+I6(2);
      end
      
      LL=[P1,P2,P3]';
      P=[P,LL]; 
  end
  PP=[p(1:2,:),s]; NN=[t(1:3,:);P];
  
  elm_bry=zeros(1,num.elem);
  
  
for i=1:num.elem
    for j=1:6
        if Nb(NN(j,i))==1
            elm_bry(i)=1;
        end
     end
end

return
end

function I=find_col(A,b)
n=length(A(1,:)); I=0;
for i=1:n
    if norm(A(:,i)-b,2)==0    
        I=[I,i];
    end    
end
return
end