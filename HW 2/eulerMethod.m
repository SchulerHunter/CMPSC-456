function x0 = eulerMethod(x0,f,j,N)
h = 1/N;
for k=1:N-1
      x0 = x0 + h * j(x0(1), x0(2));
end
end

