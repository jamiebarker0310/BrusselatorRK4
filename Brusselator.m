A = 2;
B = 6;
x0 = 0;
y0 = 1;
t0 = 0;
tend = 25;
h = 10^-1;

%Part2
[x, y, t] = brusselator(x0, y0, t0, tend, h, A, B);

%Part3 IC Inside Limit Cycle
[x2, y2, t2] = brusselator(3, 4, t0, tend, h, A, B);

%Part5 Order



tiledlayout(4,1)

%X,Y against T plot
nexttile
plot(t, x, t,y)
title('x and y as a function of t')
xlabel('t')
legend('x', 'y')

%X,Y Phase Plot
nexttile
plot(x,y, x2, y2)
legend('outside limit cycle', 'inside limit cycle')
title('phase plane of x and y')
xlabel('x')
ylabel('y')

%Order Plot
nexttile

h_vals = logspace(log10(2^-14), log10(2^-3), 12);

x_vals = zeros(11,1);
y_vals = zeros(11,1);
for i = 1:12
    [x, y, t] = brusselator(x0, y0, t0, tend, h_vals(i), A, B);
    x_vals(i) = x(end);
    y_vals(i) = y(end);
end
error_vals = (diff(x_vals).^2 + diff(y_vals).^2).^0.5;
loglog(h_vals(2:end), abs(error_vals));
m = polyfit(transpose(log(h_vals(2:end))), log(abs(error_vals)), 1)
title("gradient = "+string(m(1)))

nexttile
[x, y, t] = brusselator(x0, y0, t0, tend, 0.3, A, B);
plot(t, x)
title("Time Step = 0.3")




% error_x = abs(x2s - x(end));
% error_y = abs(y2s - y(end));
% error_norm = (error_x.^2 + error_y.^2).^0.5;
% loglog(h_vals, error_norm)
% m = polyfit(log(h_vals), log(transpose(error_norm)), 1);



function f = f(x,y,t, A, B)
    f = A - B*x + y*x^2 - x;
end

function g = g(x,y,t,A,B)
    g = B*x - y*x^2;
end

function [x, y, t] = brusselator(x0, y0, t0, t1, h, A, B)

    N = cast((t1-t0)/h, 'uint64')+1
    
    
    t = linspace(t0,t1,N);

    
    x = zeros(N,1);
    y = zeros(N,1);
    
    x(1) = x0;
    y(1) = y0;
    
    for i = 2:N
        [x(i), y(i)] = rungeKutta4(x(i-1), y(i-1), t(i-1), h, A, B);
    end
end

function [x_next, y_next] = rungeKutta4(x, y, t,h, A, B)

    k1 = f(x,y,t,A,B);
    l1 = g(x,y,t,A,B);
    
    k2 = f(x + k1*h/2, y + l1*h/2, t + h/2, A, B);
    l2 = g(x + k1*h/2, y + l1*h/2, t + h/2, A, B);
    
    k3 = f(x + k2*h/2, y + l2*h/2, t + h/2, A, B);
    l3 = g(x + k2*h/2, y + l2*h/2, t + h/2, A, B);
    
    k4 = f(x + k3*h, y+l3*h, t+h, A, B);
    l4 = g(x + k3*h, y+l3*h, t+h, A, B);
    
    x_next = x + (k1+2*k2+2*k3+k4)*h/6;
    y_next = y + (l1+2*l2+2*l3+l4)*h/6;
    
end
