%% Lecture 25 Code

% Clean workspace
clear all; close all; clc

%% Simulate and compare the steady soliton

% Space
L = 40; n = 512;
x2 = linspace(-L/2,L/2,n+1); x = x2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]';

% Time
t = 0:0.1:10;

% Initial condition
u = sech(x);
ut = fft(u);

[t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

% Obtain the POD mode
[U, S, V] = svd(usol');

% Steady soliton match
phi_xx = ifft(-(k.^2).*fft(U(:,1)));
norm = trapz(x,U(:,1).*conj(U(:,1)));
A0 = trapz(x,sech(x)'.*conj(U(:,1)))/norm;
alpha = trapz(x,phi_xx.*conj(U(:,1)))/norm;
beta = trapz(x,(U(:,1).*conj(U(:,1))).^2)/norm;

a = A0*exp(1i*0.5*alpha*t + 1i*beta*t*abs(A0)^2);

% Plot solitons
subplot(1,2,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16)
subplot(1,2,2), waterfall(x,t,abs(a*U(:,1)')), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16)

%% Compare phases of the solutions in the middle of space (index 256)

figure(2)
plot(angle(usol(:,256)),'r','Linewidth',2)
hold on
plot(angle(a*U(256,1)),'k.','Markersize',8)
xlabel('t')
ylabel('Phase')
set(gca,'FontSize',16)


%% Simulate and compare the oscillatory soliton

clear all; close all; clc

% Space
L = 40; n = 512;
x2 = linspace(-L/2,L/2,n+1); x = x2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]';

% Time
t = 0:0.1:10;

% Initial condition
u = 2*sech(x);
ut = fft(u);

[t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

% Obtain the POD modes
[U, S, V] = svd(usol');

% compute the second derivatives 
phi_1_xx = ifft( -(k.^2).*fft(U(:,1)) ); 
phi_2_xx = ifft( -(k.^2).*fft(U(:,2)) );

% compute the norms of the SVD modes 

norm1=trapz(x,U(:,1).*conj(U(:,1))); 
norm2=trapz(x,U(:,2).*conj(U(:,2)));

% compute the initial conditions 

A0_1=trapz(x,(2*sech(x).').*conj(U(:,1)))/norm1; 
A0_2=trapz(x,(2*sech(x).').*conj(U(:,2)))/norm2;

% compute inner products alpha, beta, sigma 
alpha11=0.5*trapz(x,conj(U(:,1)).*phi_1_xx)/norm1; 
alpha12=0.5*trapz(x,conj(U(:,1)).*phi_2_xx)/norm1; 
beta11_1=trapz(x,(U(:,1).*conj(U(:,1))).^2)/norm1; 
beta22_1=trapz(x,(U(:,2).*(abs(U(:,2)).^2).*conj(U(:,1))))/norm1; 
beta21_1=trapz(x,(U(:,1).*(abs(U(:,2)).^2).*conj(U(:,1))))/norm1; 
beta12_1=trapz(x,(U(:,2).*(abs(U(:,1)).^2).*conj(U(:,1))))/norm1; 
sigma12_1=trapz(x,(conj(U(:,2)).*(U(:,1).^2).*conj(U(:,1))))/norm1; 
sigma21_1=trapz(x,(conj(U(:,1)).*(U(:,2).^2).*conj(U(:,1))))/norm1;

alpha21=0.5*trapz(x,conj(U(:,2)).*phi_1_xx)/norm2; 
alpha22=0.5*trapz(x,conj(U(:,2)).*phi_2_xx)/norm2; 
beta11_2=trapz(x,(U(:,1).*(abs(U(:,1)).^2).*conj(U(:,2))))/norm2; 
beta22_2=trapz(x,(U(:,2).*(abs(U(:,2)).^2).*conj(U(:,2))))/norm2; 
beta21_2=trapz(x,(U(:,1).*(abs(U(:,2)).^2).*conj(U(:,2))))/norm2; 
beta12_2=trapz(x,(U(:,2).*(abs(U(:,1)).^2).*conj(U(:,2))))/norm2; 
sigma12_2=trapz(x,(conj(U(:,2)).*(U(:,1).^2).*conj(U(:,2))))/norm2; 
sigma21_2=trapz(x,(conj(U(:,1)).*(U(:,2).^2).*conj(U(:,2))))/norm2;

y0=[A0_1; A0_2];
[t,u_ode]=ode45(@(t,y) nls_ode_rhs(t,y,alpha11,alpha12,alpha21,alpha22, ...
beta11_1,beta22_1,beta12_1,beta21_1,beta11_2,beta22_2,beta21_2,beta12_2,...
sigma12_1,sigma21_1,sigma12_2,sigma21_2),t,y0);

    
% Plot solitons
subplot(1,2,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16,'Ylim',[0 3])
subplot(1,2,2), waterfall(x,t,abs(u_ode(:,1)*U(:,1)' + u_ode(:,2)*U(:,2)')), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16,'Ylim',[0 3])

%% Simulate a moving soliton

clear all; close all; clc

% Space
L = 40; n = 512;
x2 = linspace(-L/2,L/2,n+1); x = x2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]';

% Time
t = 0:0.1:6;

% Initial condition
u = 2*sech(x+10).*exp(1i*pi*x);
ut = fft(u);

[t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

% Solution
subplot(1,2,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u|')
set(gca,'FontSize',16)

% Solution in Fourier space
subplot(1,2,2), waterfall(k,t,abs(utsol))
xlabel('x')
ylabel('t')
zlabel('|ut|')
set(gca,'FontSize',16)

%% SVD of the moving soliton

[U, S, V] = svd(usol');

% Plot singular values
subplot(1,2,1), plot(diag(S),'k.','Markersize',20)
ylabel('\sigma_j')
set(gca,'FontSize',16,'Xlim',[0 30])

% And their logarithms
subplot(1,2,2), semilogy(diag(S),'k.','Markersize',20)
ylabel('\sigma_j')
set(gca,'FontSize',16,'Xlim',[0 30],'Ylim',[1, 1e2])

%% Find centre of moving soliton

for j = 1:length(t)
   com = x.*abs(usol(j,:)).^2; 
   com2 = abs(usol(j,:)).^2;
   c(j) = trapz(x,com)/trapz(x,com2);
end

%% Move to a moving coordinate frame

for j = 1:length(t)
   [mn,jj] = min(abs(x - c(j)));
   ns = n/2 - jj;
   ushift = [usol(j,:) usol(j,:) usol(j,:)];
   usol_shift(j,:) = ushift(1,n+1-ns:2*n-ns);
end

% Plot the result
waterfall(x,t,abs(usol_shift)), colormap([0 0 0])
xlabel('x')
ylabel('t')
zlabel('|u(x-c(t))|')
set(gca,'FontSize',16)






%% NLS Right-Hand-Side

function rhs = nls_rhs(t,ut,k)
    u = ifft(ut);
    rhs = -(1i/2)*(k.^2).*ut + 1i*fft( (abs(u).^2).*u );
end

%% NLS 2D ODE Right-Hand-Side

function rhs = nls_ode_rhs(t,y0,alpha11,alpha12,alpha21,alpha22, ...
    beta11_1,beta22_1,beta12_1,beta21_1,beta11_2,beta22_2,beta21_2,beta12_2,...
    sigma12_1,sigma21_1,sigma12_2,sigma21_2)

rhs= [ 1i*( alpha11*y0(1)+alpha12*y0(2)+(beta11_1*abs(y0(1))^2+2*beta21_1*abs(y0(2))^2)*y0(1) ... 
    +(beta12_1*abs(y0(1))^2+2*beta22_1*(abs(y0(2))^2))*y0(2) ...
    + sigma12_1*y0(1)^2*conj(y0(2)) + sigma21_1*y0(2)^2*conj(y0(1)) ); 
    1i*( alpha21*y0(1)+alpha22*y0(2)+(beta11_2*abs(y0(1))^2+2*beta21_2*abs(y0(2))^2)*y0(1) ... 
    +(beta12_2*abs(y0(1))^2+2*beta22_2*(abs(y0(2))^2))*y0(2) ...
    + sigma12_2*y0(1)^2*conj(y0(2)) + sigma21_2*y0(2)^2*conj(y0(1)) )];
end



