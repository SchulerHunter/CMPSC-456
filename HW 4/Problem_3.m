x=linspace(-1,1,100);
y=heaviside(x);
p=polyfit(x,y,25);
x1=linspace(-1,1);
p_n=polyval(p,x1);
plot(x1,p_n);