function [iter, x0] = modNewtonsMethod(x0,f,j,iMax,tol)
iter = 0;
J = j(x0(1),x0(2));
while iter < iMax
    s = linsolve(J, -f(x0(1),x0(2)));
    x0 = x0 + s;
    iter = iter + 1;
    if norm(s) < tol || norm(s) == inf
        break
    end
end
end