% Problem 2
% Using Legendre, compute coeeficients
syms phi_1(x) phi_2(x) phi_3(x) phi_4(x) f(x) p(x)
phi_1 = 1;
phi_2(x) = x;
phi_3(x) = (1/2)*(3*x^2-1);
phi_4(x) = (1/2)*(5*x^3-3*x);
f(x) = log(x+2);
a0 = integral(matlabFunction(f), -1, 1);
disp(['a0 = ', num2str(a0)])
for i = 2:4
    phi_i = eval(strcat("phi_"+num2str(i)));
    phi_i2 = phi_i*phi_i;
    q = integral(matlabFunction(f*phi_i), -1, 1);
    d = integral(matlabFunction(phi_i2), -1, 1);
    disp(['a', num2str(i),' = ', num2str(q/d)])
end
p(x) = .005*x^3 - .140*x^2 + .525*x + 1.301;
disp(['Error is proportional to ', num2str(integral(matlabFunction(f-p), -1, 1))]);