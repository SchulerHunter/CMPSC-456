clc;
syms f(x, y);
syms g(x, y);
f(x, y)=x^2-y^2; % x'
g(x, y)=x*y+x+y+1; %y'
x = zeros(1,100);
y = zeros(1,100);
x(1) = -.99999;
y(1) = .99999;
h = .01;
for i = 1:99
    % k for x; l for y
    % Solve for first variables
    k1 = h*f(x(i), y(i));
    l1 = h*g(x(i), y(i));
    % Solve for second variables
    k2 = h*f(x(i)+(k1/2), y(i)+(l1/2));
    l2 = h*g(x(i)+(k1/2), y(i)+(l1/2));
    % Solve for third variables
    k3 = h*f(x(i)+(k2/2), y(i)+(l2/2));
    l3 = h*g(x(i)+(k2/2), y(i)+(l2/2));
    % Solve for fourth variables
    k4 = h*f(x(i)+k3, y(i)+l3);
    l4 = h*g(x(i)+k3, y(i)+l3);
    % Solve for next values
    x(i+1) = x(i) + (1/3)*(.5*k1+k2+k3+.5*k4);
    y(i+1) = y(i) + (1/3)*(.5*l1+l2+l3+.5*l4);
end
plot(x, y)