function sol = contMethod(x0,f,j,N)

h = 1/N;
b = -h*f(x0(1), x0(2));
k = zeros(2,N);
sol = [0, 0];

for i = 1:N
    if i == 1
        A = j(x0(1), x0(2));
    else
        A = j(x0(1)+k(i-1, 1), x0(2)+k(i-1, 2));
    end
    A = b\A;
    k(i, 1) = A(1);
    k(i, 2) = A(2);
end


for i = 1:N
    if i == 1 || i == N
        sol(1) = sol(1) + k(i, 1);
        sol(2) = sol(2) + k(i, 2);
    else
        sol(1) = sol(1) + 2*k(i, 1);
        sol(2) = sol(2) + 2*k(i, 2);
    end
end

sol(1) = sol(1)+x0(1);
sol(2) = sol(2)+x0(2);
sol = sol./(2*N-1);
    
end

