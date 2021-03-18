%% Lecture 28 Code

% Clean workspace
clear all; close all; clc

%% Simple synthetic example 

% Generate data
dt = 0.1;
t = 0:dt:20;
x = sin(t);

% Hankel matrix
delays = 10;
xd = hankel(x(1:delays),x(delays:end));

% Apply SVD
[U, S, V] = svd(xd); % SVD of delay matrix

% Plot SVD Results
subplot(2,1,1) % Singular values
plot(diag(S),'ko','Linewidth',2)
ylabel('\sigma_j')
set(gca,'Fontsize',16,'Xlim',[0.9 delays+0.1])

subplot(2,1,2) % Right-singular vectors
plot(t(1:end-delays+1),V(:,1),'r','Linewidth',2)
hold on
plot(t(1:end-delays+1),V(:,2),'b--','Linewidth',2)
xlabel('t')
ylabel('v_j(t)')
set(gca,'Fontsize',16,'Xlim',[0 10])

%% Apply DMD to first two columns of V

X1 = V(1:end-1,1:2)';
X2 = V(2:end,1:2)'; % first two columns of V

[U2, S2, V2] = svd(X1,'econ');
Stilde = U2'*X2*V2*diag(1./diag(S2));
[eV, D] = eig(Stilde); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U2*eV;

y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions

u_modes = zeros(length(y0),length(t));
for iter = 1:length(t)
   u_modes(:,iter) = y0.*exp(omega*t(iter)); 
end
u_dmd = real(Phi*u_modes);

% print the continuous-time eigenvalues
omega

%% Simulate the Van der Pol Oscillator

% Clean workspace
clear all; close all; clc

% Simulate ODE
dt = 0.05;
t = 0:dt:200;
x0 = [2; 2];
[t, x] = ode45(@(t,x) VdP(t,x),t,x0);

% Plot solution
subplot(2,1,1) % x(t)
plot(t,x(:,1),'b','Linewidth',2)
xlabel('t')
ylabel('x(t)')
set(gca,'Fontsize',16)
axis tight

subplot(2,1,2) % y(t)
plot(t,x(:,2),'r','Linewidth',2)
xlabel('t')
ylabel('y(t)')
set(gca,'Fontsize',16)
axis tight

%% Time delay coordinates

close all

% Hankel matrix
delays = 1000;
xd = hankel(x(1:delays,1),x(delays:end,1));

% Apply SVD
[U, S, V] = svd(xd); % SVD of delay matrix

% Plot SVD Results
subplot(2,1,1) % Singular values
plot(diag(S),'ko','Linewidth',2)
ylabel('\sigma_j')
set(gca,'Fontsize',16,'Xlim',[0.9 delays+0.1])
title(['Delays = ', num2str(delays)])

subplot(2,1,2) % Right-singular vectors
plot(t(1:end-delays+1),V(:,1),'r','Linewidth',2)
hold on
plot(t(1:end-delays+1),V(:,2),'b--','Linewidth',2)
xlabel('t')
ylabel('v_j(t)')
set(gca,'Fontsize',16,'Xlim',[0 t(end-delays)])




%% VdP Right-Hand-Side

function rhs = VdP(t,x)
    rhs = [x(2); -x(1) + 10*(1 - x(1)^2)*x(2)];
end
