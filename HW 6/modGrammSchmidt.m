function Q = modGrammSchmidt(A)
%MODGRAMMSCHMIDT Uses modified grammSchmidt to create an orthonormal matrix
[rows, cols] = size(A);
Q = zeros(rows,cols);
Q(:,1) = A(:,1)/sqrt(A(:,1)'*A(:,1));
for i = 2:cols
  Q(:,i) = A(:,i);
  for j = 1:i-1
    Q(:,i) = Q(:,i) - ( Q(:,j)'*Q(:,i) )/( Q(:,j)'*Q(:,j) )*Q(:,j);
  end
  Q(:,i) = Q(:,i)/sqrt(Q(:,i)'*Q(:,i));
end
end

