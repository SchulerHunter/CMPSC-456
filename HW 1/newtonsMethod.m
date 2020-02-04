function [iter, x0] = newtonsMethod(x0,f,j,iMax,tol)
iter = 0;
while iter < iMax
    s = linsolve(j(x0(1),x0(2)), -f(x0(1),x0(2)));
    x0 = x0 + s;
    iter = iter + 1;
    if norm(s) < tol || norm(s) == inf
        break
    end
end
end