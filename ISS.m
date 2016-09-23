function [x,lambda,u,t] = ISS(L,M,f,r,x_0,T,N,m)
%ISS : 
% m = number of Control Input

lambda_0 = [-0.3;0.0];% <-0.3307,0.0492> <-0.2705,0.0736> <-0.1720,0.0063>
options = optimoptions('fsolve','MaxFunEvals',1000,'Algorithm','levenberg-marquardt');
lambda_0 = fsolve(@(lambda_0)cons(lambda_0,L,M,f,r,x_0,T,N,m),lambda_0,options); 

% define Hamiltonian function based on L, f, x, lambda, and u
H = @(x, lambda, u)L(x, u) + transpose(lambda)*f(x, u);
[x,lambda,u] = forSimIn(H,f,x_0,lambda_0,T,N,m);
t = linspace(0,T,N+1);

end

function F = cons(lambda_0,L,M,f,r,x_0,T,N,m)
% this is the equality constraint of fsolve
% this will be called in fsolve iteratively

m = 1;
l = 1;
g = 9.81;
b = 1;

H = @(x, lambda, u)L(x, u) + transpose(lambda)*f(x, u);
[x_sim,~,~] = forSimIn(H,f,x_0,lambda_0,T,N,m);
F = r(x_sim(:,end),T);%F = [x_sim(1, N+1), x_sim(2, N+1)]; % ignore the Mayer constraint for this problem TODO
end

function [x,lambda,u] = forSimIn(H,f,x_0,lambda_0,T,N,m)

dt = T/N; % step-size
n = 2; % number of states

% allocate space for x, lambda, and u
x = zeros(n,N+1);
lambda = zeros(n,N+1);
u = zeros(m,N+1);

x(:,1) = x_0;
lambda(:,1) = lambda_0;

u_max = 5;  
u(:,1) = max(-u_max,min(u_max,-lambda(2,1)/2)); 


for i = 1:N,
    y_i = [x(:,i);lambda(:,i)];
    y_next = rk4(@(y,u)dynamicAdjoint(H,f,y,u,n),y_i,u(:,i),dt);
    
    x(:,i+1) = y_next(1:n);
    lambda(:,i+1) = y_next(n+1:2*n);
    
    u(:,i+1) = max(-u_max,min(u_max,-lambda(2,i+1)/2));
end

end
function dy = dynamicAdjoint(H,f,y,u,n)
m = 1;
l = 1;
g = 9.81;
b = 1;
x = y(1:n);
lambda = y(n+1:2*n);
% determine the augmented dynamics for y ,(x lambda)
% use the given estimateGradient
dy = [f(x,u) ; -estimateGradient(@(x)H(x, lambda, u), x)];
  

end


%{
x(2) ; ...
      ((-b*x(2) + g*l*m*sin(x(1)) + u ) / (l*l*m)) ;...
      -(lambda(2)*g*cos(x(1))) / l ;...
      lambda(2)*b / (l*l*m) - lambda(1)];
 This is system-specific in this case!
%} 



function x_next = rk4(f,x,u,dt)
k1 = f(x,u);
k2 = f(x+k1*dt/2,u);
k3 = f(x+k2*dt/2,u);
k4 = f(x+k3*dt,u);
x_next = x + (k1+2*k2+2*k3+k4)*dt/6;
end