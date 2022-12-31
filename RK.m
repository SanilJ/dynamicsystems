clearvars; close all; clc;

ic = [1;2;3;4];
save('initial_conditions.mat'); 
load('initial_conditions.mat');

matObj = matfile('adjacency_matrices.mat'); 
A1 = matObj.A1; 
A2 = matObj.A2; 
A3 = matObj.A3;

L1 = diag(sum(A1,2))-A1;        
L2 = diag(sum(A2,2))-A2; 
L3 = diag(sum(A3,2))-A3; 

h = 0.5;        % step size
n = 20;         % number of steps
t0 = 0;         % initial time
y0 = ic;    % initial condition

%%% run Runge-Kutta
[t_RK1,y_RK1] = RKf(t0,y0,n,h,L1);
[t_RK2,y_RK2] = RKf(t0,y0,n,h,L2);
[t_RK3,y_RK3] = RKf(t0,y0,n,h,L3);

% plot the solution for "x" (see example 11.3)
figure,
nexttile
plot(t_RK1, y_RK1,'DisplayName','A1')
title('Runge Kutta for A1')
xlabel('timestep (t)'), ylabel('value (x)')
nexttile
plot(t_RK2, y_RK2,'DisplayName','A2')
title('Runge Kutta for A2')
xlabel('timestep (t)'), ylabel('value (x)')
nexttile
plot(t_RK3, y_RK3,'DisplayName','A3')
title('Runge Kutta for A3')
xlabel('timestep (t)'), ylabel('value (x)')


% define the function in the diff eq you want to solve: y'=f(t,y)
function ydot = f(~,y,L)
    ydot = -L*y;
end

%Runga Kutta Function from notes
function [t,y] = RKf(t0,y0,n,h,L)
t(1) = t0;      % initialize t
y(:,1) = y0;    % initialize y
for k = 1:n     % go from 1 to n to get n+1 time steps, like 0 to n in defintion of algorithm
    % incrememnt time forward
    t(k+1) = t(k)+h;
    % find y_k+1 with the Runge Kutta algorithm 
    f1 = f(t(k),y(:,k),L);
    f2 = f(t(k)+h/2,y(:,k)+h*f1/2,L);
    f3 = f(t(k)+h/2,y(:,k)+h*f2/2,L);
    f4 = f(t(k)+h,y(:,k)+h*f3,L);
    y(:,k+1) = y(:,k)+(h/6)*(f1+2*f2+2*f3+f4);
end
end
