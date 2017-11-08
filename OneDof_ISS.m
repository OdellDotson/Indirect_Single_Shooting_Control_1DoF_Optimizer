% fill inthis part
L = @(x, u) u.^2;
M = @(x,T) 0;
r = @(x,T) x; %@(x,T) fsolve(cons(lambda_0,L,M,f,r,x_0,T,N,m), x);

m = 1;
l = 1;
g = 9.81;
b = 1;

% define your dynamics
f = @(x,u,t) [x(2) ; (g/l)*sin(x(1)) - (b/(m*l*l))*x(2) + ((1/(m*l*l))*u)];
x_0 = [pi/2;0];

T = 12;
N = 80;

numControlInput = 1;
 
[x,lambda,u,t] = ISS(L,M,f,r,x_0,T,N,numControlInput);

%% Visualization

figure
dt = T/N;
step = 1;
t_previous = -step*dt;
for i = 1:step:length(t),
    theta = x(1,i);
    plot([0 -l*sin(theta)],[0 l*cos(theta)],'Linewidth',2);
    axis([-1 1 -1 1]);
    axis equal
    pause(t(i)-t_previous);
    t_previous = t(i);
    hold off;
    
end

figure
subplot(2,1,1)
plot(t,x);
xlabel('t')
ylabel('state x')
subplot(2,1,2)
plot(t,u);
xlabel('t')
ylabel('u')

