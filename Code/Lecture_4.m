%% Lecture 4 Code

% Clean workspace
clear all; close all; clc


%% Create a signal with high frequencies at the beginning and end, and low frequencies in the middle

L = 10; n = 2048;
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Create signal
S = (3*sin(2*t) + 0.5*tanh(0.5*(t-3)) + 0.2*exp(-(t-4).^2) ...
    + 1.5*sin(5*t) + 4*cos(3*(t-6).^2))/10 + (t/20).^3;
St = fft(S);

% Plot signal in both time and frequency domain
figure(1)
subplot(2,1,1) % time domain
plot(t,S,'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('S(t)')

subplot(2,1,2) % frequency domain
plot(ks,abs(fftshift(St))/max(abs(St)),'r','Linewidth',2); axis([-50 50 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('fft(S)')

%% Breaking into subdomains

figure(2)
subplot(3,2,[1 2])
plot(t,S,'k','Linewidth',2)
hold on
for t_loc = [2.5,5,7.5]
    plot([t_loc t_loc],[-1,1],'k:','Linewidth',2)
end
set(gca,'Fontsize',16)

subplot(3,2,3)
k1 = (2*pi/L)*[0:n/8-1 -n/8:-1];
k1s = fftshift(k1);
S1 = S(1:n/4);
S1t = fft(S1);
plot(k1s,abs(fftshift(S1t))/max(abs(S1t)),'r','Linewidth',2); axis([-15 15 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('fft(S)')

subplot(3,2,4)
S2 = S(n/4+1:n/2);
S2t = fft(S2);
plot(k1s,abs(fftshift(S2t))/max(abs(S2t)),'r','Linewidth',2); axis([-15 15 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('fft(S)')

subplot(3,2,5)
S3 = S(n/2+1:3*n/4);
S3t = fft(S3);
plot(k1s,abs(fftshift(S3t))/max(abs(S3t)),'r','Linewidth',2); axis([-15 15 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('fft(S)')

subplot(3,2,6)
S4 = S(3*n/4+1:end);
S4t = fft(S4);
plot(k1s,abs(fftshift(S4t))/max(abs(S4t)),'r','Linewidth',2); axis([-15 15 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('fft(S)')

%% Picture of a filter

figure(3)
plot(t,S,'k','Linewidth',2)
hold on
filter = exp(-(t-4).^4);
plot(t,filter,'m','Linewidth',3);
plot([0 4],[1.2 1.2],'m','Linewidth',3)
plot([3.2 4.8],[1.1 1.1],'m','Linewidth',3)
plot([3.2 3.2],[1.05 1.15],'m','Linewidth',3)
plot([4.8 4.8],[1.05 1.15],'m','Linewidth',3)
plot([4 4],[-1 2],'--m','Linewidth',2)
text(1.8,1.1,'\tau','Fontsize',24,'Color','m')
text(4.2,1.2,'a','Fontsize',24,'Color','m')
set(gca,'Fontsize',16), xlabel('time (t)'), ylabel('S(t)')
xticks([]), yticks([])
axis([0 10 -1 1.3])

%% Exploring different window sizes

tau = 5; % centre of window
a = 10; % window size
g = exp(-a*(t-tau).^2); % Gaussian
Sf = g.*S;
Sft = fft(Sf);

figure(4)
subplot(2,1,1) % Time domain
plot(t,Sf,'k','Linewidth',2) 
set(gca,'Fontsize',16), xlabel('Time (t)'), ylabel('S(t)*g(t-\tau)')

subplot(2,1,2) % Fourier domain 
plot(ks,abs(fftshift(Sft))/max(abs(Sft)),'r','Linewidth',2); axis([-50 50 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('FFT(S(t)*g(t-\tau))')


