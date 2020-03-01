f = @(x,y)[x^2-y^2+2*y;2*x+y^2-6];
j = @(x,y)[2*x, -2*y+2; 2, 2*y];
N = 2;
x0 = [0, 0];
contMethod(x0,f,j,N)
eulerMethod(x0,f,j,N)
x0 = [1, 1];
contMethod(x0,f,j,N)
eulerMethod(x0,f,j,N)
x0 = [3, -2];
contMethod(x0,f,j,N)
eulerMethod(x0,f,j,N)