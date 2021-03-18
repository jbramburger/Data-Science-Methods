%% Lecture 27 Code

% Clean workspace
clear all; close all; clc

%% Generate synthetic data from NLS

% Space
L = 40; n = 512;
x2 = linspace(-L/2,L/2,n+1); x = x2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]';

% Time
slices = 20;
t = linspace(0,2*pi,slices+1); dt = t(2) - t(1);

% Initial condition
u = 2*sech(x);
ut = fft(u);

% Simulating the NLS
[t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

%% Create DMD matrices

X = usol';
X1 = X(:,1:end-1);
X2 = X(:,2:end);

%% SVD of X1 and Computation of ~S

[U, Sigma, V] = svd(X1,'econ');
S = U'*X2*V*diag(1./diag(Sigma));
[eV, D] = eig(S); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U*eV;

%% Create DMD Solution

y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions

u_modes = zeros(length(y0),length(t));
for iter = 1:length(t)
   u_modes(:,iter) = y0.*exp(omega*t(iter)); 
end
u_dmd = Phi*u_modes;

%% Compare DMD and PDE Solutions

% PDE Solution
subplot(2,1,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
title('PDE Solution')
set(gca,'FontSize',16)

% DMD Solution
subplot(2,1,2), waterfall(x,t,abs(u_dmd')), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
title('DMD Solution')
set(gca,'FontSize',16)

%% Plot DMD modes

close all
waterfall(x,1:slices,abs(Phi')), colormap([0 0 0])
xlabel('x')
ylabel('modes')
zlabel('|u|')
title('DMD Modes')
set(gca,'FontSize',16)

%% Plotting Eigenvalues (omega)

% make axis lines
line = -15:15;

plot(zeros(length(line),1),line,'k','Linewidth',2) % imaginary axis
hold on
plot(line,zeros(length(line),1),'k','Linewidth',2) % real axis
plot(real(omega)*dt,imag(omega)*dt,'r.','Markersize',15)
xlabel('Re(\omega)')
ylabel('Im(\omega)')
set(gca,'FontSize',16,'Xlim',[-1.5 0.5],'Ylim',[-3 3])

%% Forecast far into the future and compute error

t = 0:dt:100;

u_modes = zeros(length(y0),length(t));
for iter = 1:length(t)
   u_modes(:,iter) = y0.*exp(omega*t(iter)); 
end
u_dmd = Phi*u_modes;

% DMD Solution
subplot(2,1,1), waterfall(x,t(1:100),abs(u_dmd(:,1:100)')), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
title('DMD Solution')
set(gca,'FontSize',16)

% Far in the future
subplot(2,1,2), waterfall(x,t(200:end),abs(u_dmd(:,200:end)')), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
title('DMD Solution')
set(gca,'FontSize',16)

%% NLS Right-Hand-Side

function rhs = nls_rhs(t,ut,k)
    u = ifft(ut);
    rhs = -(1i/2)*(k.^2).*ut + 1i*fft( (abs(u).^2).*u );
end