clc;
y = .1;
for i = 1:10
    y = solve(x*(20*x-19)-y);
    y = y(2)
end