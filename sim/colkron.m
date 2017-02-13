% COLKRON -- column-wise Kronecker product function
%-------------------------------------------------------------
% Usage: C=colkron(B,A);
%-------------------------------------------------------------
% Inputs: A -- (n x k) matrix
%         B -- (m x k) matrix
%-------------------------------------------------------------
% Output: C -- (m*n x k) matrix whose columns are the 
%              Kronecker products of the respective columns
%              of A and B.
%
%         C(:,i)=kron(B(:,i),A(:,i)); for i=1:k
%-------------------------------------------------------------
% Jim Ward
% 5/22/96
%-------------------------------------------------------------

function C=colkron(B,A)
[m,k]=size(B); [n,k]=size(A);
C=B(kron((1:m)',ones(n,1)),:).*A(kron(ones(m,1),(1:n)'),:);
return