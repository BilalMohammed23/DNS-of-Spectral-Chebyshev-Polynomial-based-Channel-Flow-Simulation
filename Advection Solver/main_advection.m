clear 
close all
clc

N = 15;
j = 0:N;
theta = j*pi/N;
x = cos(theta);

c = 1;
T = 1;

x0 = -0.3;
sigma = 0.15;
u0 = exp(-((x - x0).^2 / (sigma^2)));

u0(end) = 0;

a_t = cheb_coeff_RK4(u0, theta);    %finding a that is Chebyshev coefficient

a_t = enforce_bc_left_coeff_RK4(a_t);   %enforcing bound_cond to get correct b values

CFL = 1;
dt  = CFL * (1/N^2)/c;
nSteps = ceil(T/dt);
dt = T/nSteps;

%% Time integration (RK4 on coefficients)
for n = 1:nSteps
% Stage 1
k1 = rhs_coeff_RK4(a_t, c);

% Stage 2
a2 = a_t + 0.5*dt*k1;
a2 = enforce_bc_left_coeff_RK4(a2);
k2 = rhs_coeff_RK4(a2, c);

% Stage 3
a3 = a_t + 0.5*dt*k2;
a3 = enforce_bc_left_coeff_RK4(a3);
k3 = rhs_coeff_RK4(a3, c);

% Stage 4
a4 = a_t + dt*k3;
a4 = enforce_bc_left_coeff_RK4(a4);
k4 = rhs_coeff_RK4(a4, c);

% RK4 update
a_t = a_t + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

% Enforce BC at end of step too
a_t = enforce_bc_left_coeff_RK4(a_t);
end

u_num = cheb_eval_series_RK4(a_t, x);

figure(1)
plot(x, u0, 'o-', 'LineWidth', 1.2); hold on
plot(x, u_num, 's-', 'LineWidth', 1.8);
grid on
legend('Time = 0', 'Time = T')

xlim([-1,1])
ylim([-0.2,1.2])
xlabel('x'); ylabel('u')
title(sprintf('Advection on [-1,1] with inflow BC u(-1,t)=0 (N=%d, T=%.2f)', N, T))