f = @(x, xp) -1*(x+xp^2);
fx = @(x, xp) -1;
fxp = @(x, xp) -2*xp;

x0 = -25;
xp0 = 6;

h = .1;
a = 0;
b = 30;

x = ab3(f, fx, fxp, x0, xp0, h, a, b);
plot(x(1, :),x(2, :));

f = @(x, xp) -sin(xp);
fx = @(x, xp) 0;
fxp = @(x, xp) cos(xp);

x0 = 1;
x1 = -1;
b = 1;
iter = 4;

x = shooting(f, fx, fxp, x0, x1, h, a, b, iter);
disp(x');

function x = ab3(f, fx, fxp, x0, xp0, h, a, b)
    N = (b-a)/h;
    x = zeros(2, N+1);
    k = zeros(4, 2);
    kp = zeros(4, 2);
    
    x(1, 1) = x0;
    x(2, 1) = xp0;
    u1 = 0;
    u2 = 1;
    
    R1 = @(x,k) .5*k+x;
    R2 = @(mat,i) (1/6)*(mat(1,i)+2*mat(2,i)+2*mat(3,i)+mat(4,i));
    % Run Runge-Kutta method for starting data
    for i = 1:2
        k(1, 1) = h*x(2, i);
        k(1, 2) = h*f(x(1,i), x(2,i));
        for j = 2:3
            k(j, 1) = h*R1(x(2,i),k(j-1,2));
            k(j, 2) = h*f(R1(x(1,i),k(j-1,1)),R1(x(2,i),k(j-1,2)));
        end
        k(4, 1) = h*(x(2,i)+k(3,2));
        k(4, 2) = h*f(x(1,i)+k(3,1),x(2,i)+k(3,2));
        
        x(1, i+1) = x(1,i) + R2(k,1);
        x(2, i+1) = x(2,i) + R2(k,2);
        
        fxTemp = fx(x(1,i),x(2,i));
        fxpTemp = fxp(x(1,i),x(2,i));
        
        kp(1, 1) = h*u2;
        kp(1, 2) = h*u1*fxTemp+u2*fxpTemp;
        for j = 2:3
            kp(j, 1) = h*R1(u2,kp(j-1,2));
            kp(j, 2) = h*R1(u1,kp(j-1,1))*fxTemp+R1(u2,kp(j-1,2))*fxpTemp;
        end
        kp(4, 1) = h*(u2 + kp(3,2));
        kp(4, 2) = h*(u1+kp(3,1))*fxTemp+(u2+kp(3,2));
        
        u1 = u1 + R2(kp,1);
        u2 = u2 + R2(kp,2);
    end
    
    for i = 3:N
        x(1, i+1) = x(1,i)+(h/12)*(23*x(2,i)-16*x(2,i-1)+5*x(2,i-2));
        x(2, i+1) = x(2,i)+(h/12)*(23*f(x(1,i+1),x(2,i))-16*f(x(1,i),x(2,i-1))+5*f(x(1,i-1),x(2,i-2)));
    end
end

function x = shooting(f, fx, fxp, x0, x1, h, a, b, max)
    N = (b-a)/h;
    iter = 1;
    TK = (x1-x0)/(b-a);
    tol = 10^-5;
    
    x = zeros(1, max);
    tempX = zeros(2, N+1);
    k = zeros(4, 2);
    kp = zeros(4, 2);
    
    R1 = @(x,k) .5*k+x;
    R2 = @(mat,i) (1/6)*(mat(1,i)+2*mat(2,i)+2*mat(3,i)+mat(4,i));
    while(iter <= max)
        tempX(1, 1) = x0;
        tempX(2, 1) = TK;
        u1 = 0;
        u2 = 1;
        for i = 1:N
            k(1, 1) = h*tempX(2, i);
            k(1, 2) = h*f(tempX(1,i), tempX(2,i));
            for j = 2:3
                k(j, 1) = h*R1(tempX(2,i),k(j-1,2));
                k(j, 2) = h*f(R1(tempX(1,i),k(j-1,1)),R1(tempX(2,i),k(j-1,2)));
            end
            k(4, 1) = h*(tempX(2,i)+k(3,2));
            k(4, 2) = h*f(tempX(1,i)+k(3,1),tempX(2,i)+k(3,2));

            tempX(1, i+1) = tempX(1,i) + R2(k,1);
            tempX(2, i+1) = tempX(2,i) + R2(k,2);

            fxTemp = fx(tempX(1,i),tempX(2,i));
            fxpTemp = fxp(tempX(1,i),tempX(2,i));

            kp(1, 1) = h*u2;
            kp(1, 2) = h*u1*fxTemp+u2*fxpTemp;
            for j = 2:3
                kp(j, 1) = h*R1(u2,kp(j-1,2));
                kp(j, 2) = h*R1(u1,kp(j-1,1))*fxTemp+R1(u2,kp(j-1,2))*fxpTemp;
            end
            kp(4, 1) = h*(u2 + kp(3,2));
            kp(4, 2) = h*(u1+kp(3,1))*fxTemp+(u2+kp(3,2));

            u1 = u1 + R2(kp,1);
            u2 = u2 + R2(kp,2);
            
            if(abs(tempX(1, N+1) - x1) <= tol)
                tempX(1:2, N+1) = tempX(1:2, i+1);
                break
            end
        end
        x(1, iter) = tempX(1, N+1);
        TK = TK - (tempX(1, N+1) - x1) / u1;
        iter = iter + 1;
    end
end