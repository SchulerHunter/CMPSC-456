function [vec,value]=invPowerMethod(start,A,toler)
%
%Inverse power method for computing eigenvalues
%
dd = 1;
x = start;
n = 10;
while dd> toler
    y = A\x; % Same as inv(A)*x
    n = max(abs(y))
    x = y/n;
    dd = norm(A\x-n*x);
end
vec=x;
value=n;