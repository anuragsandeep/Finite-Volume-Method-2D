%% Tridiagonal matrix algorithm
%% Anurag Sandeep K. 
% A=co-efficient matrix
% B=solution vector
% N=number of equations/nodes

function [B] = Tridiagonal(N,A,B)

 A(1,3)=-A(1,3)/A(1,2);
 B(1)=B(1)/A(1,2);
 for i=2:N
 A(i,3)=-A(i,3)/(A(i,2)+A(i,1)*A(i-1,3));
 B(i)=(B(i)-A(i,1)*B(i-1))/(A(i,2)+A(i,1)*A(i-1,3));
 end
 for i=N-1:-1:1
 B(i)=A(i,3)*B(i+1)+B(i);
 
end
