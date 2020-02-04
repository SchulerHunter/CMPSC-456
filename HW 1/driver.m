f = @(x,y)[x^2-y^2+2*y; 2*x+y^2-6];
j = @(x,y)[2*x, -2*y+2; 2, 2*y];
x0 = [-5, -4];
iMax = 100;
tol = 1;
newtonsMethod(x0,f,j,iMax,tol)
modNewtonsMethod(x0,f,j,iMax,tol)