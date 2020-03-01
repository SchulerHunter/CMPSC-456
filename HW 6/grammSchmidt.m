function [Q, R] = grammSchmidt(A)
%GRAMMSCHMIDT Create a Q for QR factorization
%   A function which uses the Gramm-Schmidt process to create an orthgonal
%   Q matrix for QR factorization
%   A is an n*m column matrix
%   Q is an orthogonal n*m matrix
%   R is is an m*n upper triangular matrix
[rows, cols] = size(A);
R = zeros(cols);
R(1,1)=norm(A(:,1));Q(:,1)=A(:,1)/R(1,1);
for k=2:cols
    R(1:k-1,k)=Q(:,1:k-1)'*A(:,k);
    Q(:,k)=A(:,k)-Q(:,1:k-1)*R(1:k-1,k);
    R(k,k)=norm(Q(:,k));
    Q(:,k)=Q(:,k)/R(k,k);
end
end
