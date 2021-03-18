%% Lecture 1 Code

% Clean workspace
clear all; close all; clc


%% Play audio sample and plot in time domain
load handel % a clip of Handel's "Messiah" is included with MATLAB

v = y'/2;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');

soundsc(v,Fs) % Play the audio

%% Fourier series for continuous vs discontinous functions

% Continuous function: f(x) = x^2
x = linspace(-pi,pi,10001);
y = x.^2;

% Plot x^2 function
figure(1)
plot(x,y,'b','Linewidth',4)
set(gca, 'Fontsize', 16)
xlabel('x')
ylabel('y')

% Plot Nth partial sum of Fourier series
N = 10;
yN = pi^2/3*ones(size(x)); %0th term
for k = 1:N %1 to Nth terms
    yN = yN + (-1)^k*4/(k^2)*cos(k*x);
end

hold on
plot(x,yN,'r','Linewidth',2)
xlim([-pi,pi])

% Discontinuous function: unit step function
y = (x > 0);

figure(2)
plot(x(1:5001),y(1:5001),'b','Linewidth',2)
hold on
plot([0 0],[0 1],'--b','Linewidth',4)
plot(x(5002:10000),y(5002:10000),'b','Linewidth',2)
axis([-pi pi -.2 1.2])

% Plot Nth partial sum of Fourier series
N = 200;
yN = 0.5*ones(size(x));
for k = 1:N
   yN = yN + 2/((2*k-1)*pi)*sin((2*k-1)*x); 
end

plot(x,yN,'r','Linewidth',1)
set(gca, 'Fontsize', 16)
xlabel('x')
ylabel('y')

%% Aliasing

N = 8;
x2 = linspace(-pi,pi,N+1); % domain discretization
x = x2(1:N); % consider only the first N points (periodicity)
xplot = linspace(-pi,pi,10001); %domain for plotting

% sin((N-1)*x) function plot
y1 = sin((N-1)*x);
y1plot = sin((N-1)*xplot);
plot(xplot,y1plot,'b','Linewidth',2)
xlabel('x')
set(gca, 'Fontsize', 16)

% Add in sin(-x)
hold on
y2plot = sin(-xplot);
plot(xplot,y2plot,'r','Linewidth',2)

% Add in sin((2*N-1)*x)
y3plot = sin((2*N-1)*xplot);
plot(xplot,y3plot,'g','Linewidth',2)
xlim([-pi,pi])

% Add in discrete points
plot(x,y1,'ok','Markersize',10)

%% Calculate the fft of a Gaussian

L = 20;
N = 128; % N = 2^7

x2 = linspace(-L/2,L/2,N+1); % domain discretization
x = x2(1:N); % consider only the first N points (periodicity)

u = exp(-x.^2); % Gaussian function

subplot(3,1,1)
plot(x,u,'Linewidth',2)
set(gca, 'Fontsize', 16)
xlabel('x')
ylabel('u')

% Plot fft of u
ut = fft(u); % fft of u
k = -N/2:N/2-1; %frequency domain
subplot(3,1,2)
plot(k,abs(ut),'r','Linewidth',2)
set(gca, 'Fontsize', 16)
xlabel('k')
ylabel('ut')
xlim([-64,63])
% the frequencies and ut don't match up because ut has frequencies in the
% order k = 0,1,...,N/2-1,-N/2,...,-1

% Plot shifted fft of u
utshift = fftshift(ut); % shift Fourier transform ut
subplot(3,1,3)
plot(k,abs(utshift),'r','Linewidth',2)
set(gca, 'Fontsize', 16)
xlabel('k')
ylabel('ut')
xlim([-64,63])

%% Differentiation: spectral method vs finite differences

L = 20; % length of the computational domain: [-L/2,L/2]
n = 128; % number of Fourier modes: n = 2^7

x2 = linspace(-L/2,L/2,n+1); %domain discretization
x = x2(1:n); % only the first n points (periodicity)
dx = x(2) - x(1); % finite difference step size
u = sech(x); % function to take the derivative of
ut = fft(u); % fft of the function 
k = (2*pi/L)*[0:(n/2 - 1) (-n/2):-1]; % rescaled frequencies

% fft derivatives
dutdk = 1i*k.*ut; % first derivative
d2utdk2 = -k.^2.*ut; % second derivative
dudx = real(ifft(dutdk)); %inverse transform - round off imaginary part due to error
d2udx2 = real(ifft(d2utdk2));

dudx_exact = -sech(x).*tanh(x); % analytic first derivative
d2udx2_exact = sech(x)-2*sech(x).^3; % analytic second derivative

% Finite difference second derivative

% Second order accuracy
ufd2(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx); 
for j=2:n-1
    ufd2(j)=(u(j+1)-u(j-1))/(2*dx); 
end
ufd2(n)=(3*u(n)-4*u(n-1)+u(n-2))/(2*dx); 

% Fourth order accuracy 
ufd4(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx); 
ufd4(2)=(-3*u(2)+4*u(3)-u(4))/(2*dx); 
for j=3:n-2
    ufd4(j)=(-u(j+2)+8*u(j+1)-8*u(j-1)+u(j-2))/(12*dx); 
end
ufd4(n-1)=(3*u(n-1)-4*u(n-2)+u(n-3))/(2*dx);
ufd4(n)=(3*u(n)-4*u(n-1)+u(n-2))/(2*dx);

% Spectral derivatives
figure(1)
plot(x,u,'r',x,dudx,'g',x,dudx_exact,'go',x,d2udx2,'b',x,d2udx2_exact,'bo','Linewidth',2,'Markersize',10)
legend('u(x)',"Spectral u'(x)","Exact u'(x)","Spectral u''(x)","Exact u''(x)",'Fontsize',16)
xlabel('x')
set(gca,'Fontsize',16)

% Compare first derivatives
figure(2)
plot(x,dudx_exact,'bs-',x,dudx,'gv',x,ufd2,'ro',x,ufd4,'c*','Linewidth',2,'Markersize',10)
legend('Exact','fft','2nd order FD','4th order FD','Fontsize',16)
xlabel('x')
title("u'(x)",'Fontsize',16)
set(gca,'Fontsize',16)

% Zoom in on first derivatives
figure(3)
subplot(3,1,1), plot(x,dudx_exact,'bs-',x,dudx,'gv',x,ufd2,'ro',x,ufd4,'c*')
set(gca,'Fontsize',12)
axis([-1.15 -0.75 0.47 0.5])
set(gca,'Fontsize',16)
subplot(3,1,2); plot(x,dudx_exact,'bs-',x,dudx,'gv',x,ufd2,'ro',x,ufd4,'c*') 
axis([-0.9376 -0.9374 0.49848 0.49850])
set(gca,'Fontsize',16)
subplot(3,1,3), plot(x,dudx_exact,'bs-',x,dudx,'gv',x,ufd2,'ro',x,ufd4,'c*') 
axis([-0.9376 -0.9374 0.498487 0.498488])
set(gca,'Fontsize',16)

% ------------------------------------------------------------------------
% Function with a jump discontinuity at the boundaries
u = tanh(x); % function to take a derivative of
ut = fft(u); % fft of the function

% FFT calculation of derivatives
dutdk = 1i*k.*ut; % first derivative
dudx = real(ifft(dutdk)); %inverse transform - round off imaginary part due to error
dudx_exact = sech(x).^2; % analytic first derivative

% 4th-order accurate finite difference

ufd2(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx); 
ufd2(2)=(-3*u(2)+4*u(3)-u(4))/(2*dx); 
for j=3:n-2
    ufd2(j)=(-u(j+2)+8*u(j+1)-8*u(j-1)+u(j-2))/(12*dx); 
end
ufd2(n-1)=(3*u(n-1)-4*u(n-2)+u(n-3))/(2*dx);
ufd2(n)=(3*u(n)-4*u(n-1)+u(n-2))/(2*dx);

figure(4)
plot(x,u,'b',x,dudx_exact,'ko',x,dudx,'r',x,ufd2,'g.','Linewidth',2,'Markersize',6) 
axis([-10 10 -1.5 1.5])
legend('u(x)',"u'(x)",'fft','4th order FD','Fontsize',16,'Location','Best')
xlabel('x')
title("u(x) = tanh(x)",'Fontsize',16)
set(gca,'Fontsize',16)