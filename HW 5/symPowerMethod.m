function [vec,value]=symPowerMethod(start,A,toler)
%
%Power method for computing eigenvalues of symmetric matrices
%
dd = 1;
x = start/norm(start);
n = 10;
while dd> toler
    y = A*x;
    n = y.'*x
    dd = norm((y/norm(y))-x);
    x = y/norm(y);
end
vec=x;
value=n;