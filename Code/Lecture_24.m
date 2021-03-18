%% Lecture 24 Code

% Clean workspace
clear all; close all; clc

%% Initializations for NLS

% Space
L = 40; n = 512;
x2 = linspace(-L/2,L/2,n+1); x = x2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]';

% Time
t = 0:0.1:10;

% Initial condition
u = sech(x);
ut = fft(u);

%% Simulating the NLS

[t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

% Solution
subplot(2,1,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16)

% Solution in Fourier space
subplot(2,1,2), waterfall(k,t,abs(utsol))
xlabel('x')
ylabel('t')
zlabel('|ut|')
set(gca,'FontSize',16)

%% POD modes of solution

% SVD
[U, S, V] = svd(usol);

% Plot singular values
subplot(1,2,1), plot(diag(S),'k.','Markersize',20)
ylabel('\sigma_j')
set(gca,'FontSize',16,'Xlim',[0 30])

% And their logarithms
subplot(1,2,2), semilogy(diag(S),'k.','Markersize',20)
ylabel('\sigma_j')
set(gca,'FontSize',16,'Xlim',[0 30],'Ylim',[1e-7, 1e3])

%% Different Initial condition

% Initial condition
u = 2*sech(x);
ut = fft(u);

[t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

% Solution
subplot(2,1,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16)

% Solution in Fourier space
subplot(2,1,2), waterfall(k,t,abs(utsol))
xlabel('x')
ylabel('t')
zlabel('|ut|')
set(gca,'FontSize',16)

%% POD modes of the solution

% SVD
[U, S, V] = svd(usol);

% Plot singular values
subplot(1,2,1), plot(diag(S),'k.','Markersize',20)
ylabel('\sigma_j')
set(gca,'FontSize',16,'Xlim',[0 30])

% And their logarithms
subplot(1,2,2), semilogy(diag(S),'k.','Markersize',20)
ylabel('\sigma_j')
set(gca,'FontSize',16,'Xlim',[0 30],'Ylim',[1e-7, 1e3])


%% NLS Right-Hand-Side

function rhs = nls_rhs(t,ut,k)
    u = ifft(ut);
    rhs = -(1i/2)*(k.^2).*ut + 1i*fft( (abs(u).^2).*u );
end



