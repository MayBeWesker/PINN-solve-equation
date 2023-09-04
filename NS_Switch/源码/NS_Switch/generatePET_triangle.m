function [p,e,t]=generatePET_triangle(left,right,bottom,top,h_partition)

basis_type=1;

N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);

N1_basis=N1_partition;
N2_basis=N2_partition;

[p,t]=generate_M_T_triangle(left,right,bottom,top,h_partition,basis_type);

[boundary_nodes,e]=generate_boundary_nodes_edges(N1_basis,N2_basis,N1_partition,N2_partition);

%[boundary_nodes,e]=generate_boundary_nodes_edges_Stokes(N1_basis,N2_basis,N1_partition,N2_partition);

return
end

function [M,T]=generate_M_T_triangle(left,right,bottom,top,h_partition,basis_type)

h=h_partition;

if basis_type==1

   N1=(right-left)/h(1);
   N2=(top-bottom)/h(2);
   tnp=(N1+1)*(N2+1);
   M=zeros(2,tnp);
   T=zeros(3,2*N1*N2);
   Q=zeros(N1+1,N2+1);

   for j=1:tnp
      if mod(j,N2+1)==0
         M(1,j)=left+(j/(N2+1)-1)*h(1);
         M(2,j)=top;
      else
         M(1,j)=left+fix(j/(N2+1))*h(1);
         M(2,j)=bottom+(mod(j,N2+1)-1)*h(2);
      end
   end

   for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,2*n-1)=Q(column,row);
      T(2,2*n-1)=Q(column+1,row);
      T(3,2*n-1)=Q(column,row+1);  
  
      T(1,2*n)=Q(column,row+1);
      T(2,2*n)=Q(column+1,row);
      T(3,2*n)=Q(column+1,row+1);  
    
   end

elseif basis_type==2

   N1=(right-left)/h(1);
   N2=(top-bottom)/h(2);
   dh=h/2;
   dN1=N1*2;
   dN2=N2*2;
   tnp=(dN1+1)*(dN2+1);
   M=zeros(2,tnp);
   T=zeros(3,2*N1*N2);
   Q=zeros(dN1+1,dN2+1);

   for j=1:tnp
      if mod(j,dN2+1)==0
         M(1,j)=left+(j/(dN2+1)-1)*dh(1);
         M(2,j)=top;
      else
         M(1,j)=left+fix(j/(dN2+1))*dh(1);
         M(2,j)=bottom+(mod(j,dN2+1)-1)*dh(2);
      end
   end

   for i=1:dN1+1
      for j=1:dN2+1
         Q(i,j)=(i-1)*(dN2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,2*n-1)=Q(2*column-1,2*row-1);
      T(2,2*n-1)=Q(2*column+1,2*row-1); 
      T(3,2*n-1)=Q(2*column-1,2*row+1);
      T(4,2*n-1)=Q(2*column,2*row-1);
      T(5,2*n-1)=Q(2*column,2*row);
      T(6,2*n-1)=Q(2*column-1,2*row);


      T(1,2*n)=Q(2*column-1,2*row+1);
      T(2,2*n)=Q(2*column+1,2*row-1);
      T(3,2*n)=Q(2*column+1,2*row+1);
      T(4,2*n)=Q(2*column,2*row);
      T(5,2*n)=Q(2*column+1,2*row);
      T(6,2*n)=Q(2*column,2*row+1); 

   end
end
return
end

function [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges(N1,N2,N1_partition,N2_partition)

nbn=2*(N1+N2);
boundary_nodes=zeros(2,nbn);

%The following boundary condition may change for different problems.
%All Dirichlet boundary nodes.
boundary_nodes(1,:)=-1;

%The index in the following is associated with the index in "generate_M_T_triangle.m".
%bottom boundary nodes.
for k=1:N1
    boundary_nodes(2,k)=(k-1)*(N2+1)+1;
end

%right boundary nodes.
for k=N1+1:N1+N2
    boundary_nodes(2,k)=N1*(N2+1)+k-N1;
end

%top boundary nodes.
for k=N1+N2+1:2*N1+N2
    boundary_nodes(2,k)=(2*N1+N2+2-k)*(N2+1);
end

%left boundary nodes.
for k=2*N1+N2+1:nbn
    boundary_nodes(2,k)=2*N1+2*N2+2-k;
end





%Information matrix for boundary edges. It uses the index of partition, not the index of FE.

nbe=2*(N1_partition+N2_partition);
boundary_edges=zeros(4,nbe);


%The following boundary condition may change for different problems.
%All Dirichlet boundary edges.
boundary_edges(1,:)=-1;

%The index in the following is associated with the index in "generate_M_T_triangle.m".
%bottom boundary edges.
for k=1:N1_partition
    boundary_edges(2,k)=(k-1)*2*N2_partition+1;
    boundary_edges(3,k)=(k-1)*(N2_partition+1)+1;
    boundary_edges(4,k)=k*(N2_partition+1)+1;
end

%right boundary edges.
for k=N1_partition+1:N1_partition+N2_partition
    boundary_edges(2,k)=(N1_partition-1)*2*N2_partition+2*(k-N1_partition);
    boundary_edges(3,k)=N1_partition*(N2_partition+1)+k-N1_partition;
    boundary_edges(4,k)=N1_partition*(N2_partition+1)+k-N1_partition+1;
end

%top boundary edges.
for k=N1_partition+N2_partition+1:2*N1_partition+N2_partition
    boundary_edges(2,k)=(2*N1_partition+N2_partition+1-k)*2*N2_partition;
    boundary_edges(3,k)=(2*N1_partition+N2_partition+2-k)*(N2_partition+1);
    boundary_edges(4,k)=(2*N1_partition+N2_partition+1-k)*(N2_partition+1);
end

%left boundary edges.
for k=2*N1_partition+N2_partition+1:nbe
    boundary_edges(2,k)=2*(2*N1_partition+2*N2_partition+1-k)-1;
    boundary_edges(3,k)=2*N1_partition+2*N2_partition+2-k;
    boundary_edges(4,k)=2*N1_partition+2*N2_partition+1-k;
end
return
end

function [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges_Stokes(N1_basis,N2_basis,N1_partition,N2_partition)
%Xiaoming He, 07/08/2009.
%Generate the information matrices for boundary nodes and boundary edges for Stokes equations.
%This function needs to be modified for different boundary conditons.
%Here, this function is for pure Dirichlet boundary condition.
%This function needs to be modified if we change the index of nodes and elements in "generate_M_T_triangle.m".
%We will use "FE" to replace "finite element" in the comments.
%N1_basis,N2_basis:The N1 and N2 for the FE basis functions,not the partition.
%N1_partition,N2_partition: The N1 and N2 for the partition, not the FE basis functions.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.


%boundary_nodes(1,k): specifiy the type of the kth boundary node for the normal direction(or u1).
%boundary_nodes(1,k)=-1: Dirichlet boundary node in the normal direction(or u1);
%boundary_nodes(1,k)=-2: Stress boundary node in the normal direction(or u1);
%boundary_nodes(1,k)=-3: Robin boundary node in the normal direction(or u1). 
%boundary_nodes(2,k): specifiy the type of the kth boundary node for the tangential direction(or u2).
%boundary_nodes(2,k)=-1: Dirichlet boundary node in the tangential direction(or u2);
%boundary_nodes(2,k)=-2: Stress boundary node in the tangential direction(or u2);
%boundary_nodes(2,k)=-3: Robin boundary node in the tangential direction(or u2).
%The intersection node between Dirichlet boundary and other boundaries is a Dirichlet boundary node.
%boundary_nodes(3,k): global index of the kth boundary node among all nodes of FE. 
%                     That is, the index of FE is used here.
%boundary_nodes(4:5,k): the unit outer normal vector at the kth boundary node.
%boundary_nodes(6:7,k): the unit tangential vector at the kth boundary node. 
%                       It's counterclockwise to go from the normal vector to the tangential vector.

%boundary_edges(1,k): specifiy the type of the kth boundary edge in normal direction.
%boundary_edges(1,k)=-1: Dirichlet boundary edge in normal direction;
%boundary_edges(1,k)=-2: Stress boundary edge in normal direction;
%boundary_edges(1,k)=-3: Robin boundary edge in normal direction.
%boundary_edges(2,k): specifiy the type of the kth boundary edge in tangential direction.
%boundary_edges(2,k)=-1: Dirichlet boundary edge in tangential direction;
%boundary_edges(2,k)=-2: Stress boundary edge in tangential direction;
%boundary_edges(2,k)=-3: Robin boundary edge in tangential direction.
%boundary_edges(3,k): index of the element which contains the kth boundary edge.
%boundary_edges(4:5,k): indices of the two end points of the kth boundary edge among all grid points, not the nodes of FE.
%                       That is, the index of partition is used here.
%boundary_edges(6:7,k): the unit outer normal vector at the kth boundary edge.
%第k边界边上的单位外法向量。
%boundary_edges(8:9,k): the unit tangential vector at the kth boundary edge. 
%第k个边界边处的单位切向量。
%                       It's counterclockwise to go from the normal vector to the tangential vector.
%More explanation is in my "Notes for tool box of standard triangular FE" section 2-3-7.

%nbn: the total number of all the boundary nodes of FE.
%nbe: the total number of all the boundary edges.


%% Information matrix for boundary nodes. It uses the index of FE, not the index of partition.

nbn=2*(N1_basis+N2_basis);
boundary_nodes=zeros(3,nbn);

%The following boundary condition may change for different problems.
%All boundary nodes are Dirichlet boundary nodes for both normal and tangential directions(or u1 and u2).
boundary_nodes(1,:)=-1;
boundary_nodes(2,:)=-1;

%The index in the following is associated with the index in "generate_M_T_triangle.m".
%bottom boundary nodes.
for k=1:N1_basis
    boundary_nodes(3,k)=(k-1)*(N2_basis+1)+1;
    boundary_nodes(4:5,k)=[0 -1];
end

%right boundary nodes.
for k=N1_basis+1:N1_basis+N2_basis
    boundary_nodes(3,k)=N1_basis*(N2_basis+1)+k-N1_basis;
    boundary_nodes(4:5,k)=[1 0];
end

%top boundary nodes.
for k=N1_basis+N2_basis+1:2*N1_basis+N2_basis
    boundary_nodes(3,k)=(2*N1_basis+N2_basis+2-k)*(N2_basis+1);
    boundary_nodes(4:5,k)=[0 1];
end

%left boundary nodes.
for k=2*N1_basis+N2_basis+1:nbn
    boundary_nodes(3,k)=2*N1_basis+2*N2_basis+2-k;
    boundary_nodes(4:5,k)=[-1 0];
end

%Correct the normal direction at the four corners.
%left-bottom corner.
boundary_nodes(4:5,1)=[-1 -1]/sqrt(2);
%right-bottom corner.
boundary_nodes(4:5,N1_basis+1)=[1 -1]/sqrt(2);
%right-top corner.
boundary_nodes(4:5,N1_basis+N2_basis+1)=[1 1]/sqrt(2);
%left-top corner.
boundary_nodes(4:5,2*N1_basis+N2_basis+1)=[-1 1]/sqrt(2);

%It's counterclockwise to go from the normal vector n to the tangential vector \tau.
%Hence \tau_1=-n_2, \tau_2=n_1.
boundary_nodes(6,:)=-boundary_nodes(5,:);
boundary_nodes(7,:)=boundary_nodes(4,:);



%% Information matrix for boundary edges. It uses the index of partition, not the index of FE.

nbe=2*(N1_partition+N2_partition);
boundary_edges=zeros(5,nbe);


%The following boundary condition may change for different problems.
%All boundary edges are Dirichlet boundary edges for both normal and tangential directions.
boundary_edges(1,:)=-1;
boundary_edges(2,:)=-1;

%The index in the following is associated with the index in "generate_M_T_triangle.m".
%bottom boundary edges.
for k=1:N1_partition
    boundary_edges(3,k)=(k-1)*2*N2_partition+1;
    boundary_edges(4,k)=(k-1)*(N2_partition+1)+1;
    boundary_edges(5,k)=k*(N2_partition+1)+1;
    boundary_edges(6:7,k)=[0 -1];
end

%right boundary edges.
for k=N1_partition+1:N1_partition+N2_partition
    boundary_edges(3,k)=(N1_partition-1)*2*N2_partition+2*(k-N1_partition);
    boundary_edges(4,k)=N1_partition*(N2_partition+1)+k-N1_partition;
    boundary_edges(5,k)=N1_partition*(N2_partition+1)+k-N1_partition+1;
    boundary_edges(6:7,k)=[1 0];
end

%top boundary edges.
for k=N1_partition+N2_partition+1:2*N1_partition+N2_partition
    boundary_edges(3,k)=(2*N1_partition+N2_partition+1-k)*2*N2_partition;
    boundary_edges(4,k)=(2*N1_partition+N2_partition+2-k)*(N2_partition+1);
    boundary_edges(5,k)=(2*N1_partition+N2_partition+1-k)*(N2_partition+1);
    boundary_edges(6:7,k)=[0 1];
end

%left boundary edges.
for k=2*N1_partition+N2_partition+1:nbe
    boundary_edges(3,k)=2*(2*N1_partition+2*N2_partition+1-k)-1;
    boundary_edges(4,k)=2*N1_partition+2*N2_partition+2-k;
    boundary_edges(5,k)=2*N1_partition+2*N2_partition+1-k;
    boundary_edges(6:7,k)=[-1 0];
end


%It's counterclockwise to go from the normal vector n to the tangential vector \tau.
%Hence \tau_1=-n_2, \tau_2=n_1.
boundary_edges(8,:)=-boundary_edges(7,:);
boundary_edges(9,:)=boundary_edges(6,:);
return
end




