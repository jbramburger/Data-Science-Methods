%% Lecture 5 Code

% Clean workspace
clear all; close all; clc

%% Plotting the Haar wavelet and its Fourier transform

L = 8; n = 256;
t2 = linspace(-L/2,L/2,n+1); t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Create Haar wavelet
H = zeros(length(t),1);
for j = 1:length(t)
   if t(j) >= 0 && t(j) < 0.5
       H(j) = 1;
   elseif t(j) >= 0.5 && t(j) < 1
       H(j) = -1;
   end
end

% Fourier transform of Haar wavelet
Ht = fft(H);

% Plot Haar wavelet in both time and frequency domain
figure(1)
subplot(1,2,1) % time domain
plot(t,H,'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('\psi(t)')
axis([-1 2 -1.2 1.2])

subplot(1,2,2) % frequency domain
plot(ks,abs(fftshift(Ht))/max(abs(Ht)),'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('frequency (k)'); ylabel('FFT(\psi)')
axis([-100 100 0 1.05])


%% Plotting scaled and shifted Haar wavelets

% Create the scaled and shifted Haar wavelets
H1 = Haar(t,0.5,0);
H1t = fft(H1);
H2 = Haar(t,2,0);
H2t = fft(H2);
H3 = Haar(t,2,-0.5);
H3t = fft(H3);

figure(1)
subplot(3,2,1) % time domain
plot(t,H1,'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('\psi_{1/2,0}(t)')
axis([-1 3 -2 2])

subplot(3,2,2) % frequency domain
plot(ks,abs(fftshift(H1t))/max(abs(H1t)),'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('frequency (k)'); ylabel('FFT(\psi_{1/2,0})')
axis([-100 100 0 1.05])

subplot(3,2,3) % time domain
plot(t,H2,'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('\psi_{2,0}(t)')
axis([-1 3 -2 2])

subplot(3,2,4) % frequency domain
plot(ks,abs(fftshift(H2t))/max(abs(H2t)),'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('frequency (k)'); ylabel('FFT(\psi_{2,0})')
axis([-100 100 0 1.05])

subplot(3,2,5) % time domain
plot(t,H3,'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('\psi_{2,-1/2}(t)')
axis([-1 3 -2 2])

subplot(3,2,6) % frequency domain
plot(ks,abs(fftshift(H3t))/max(abs(H3t)),'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('frequency (k)'); ylabel('FFT(\psi_{2,-1/2})')
axis([-100 100 0 1.05])

%% Creating new wavelets with convolutions

% Create the Gaussian and fft
phi = exp(-t.^2)';
phit = fft(phi);

% Multiply transforms and apply ifft
convt = Ht.*phit;
conv = ifft(convt);

figure(3)
plot(t,ifftshift(conv),'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('\psi\ast\phi(t)')
axis([-4 4 -7 7])

%% Plotting the Mexican hat wavelet

L = 50; n = 1024;
t2 = linspace(-L/2,L/2,n+1); t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

mex = (1 - t.^2).*exp(-0.5*t.^2);
mext = fft(mex);

% Plot Mexican hat wavelet in both time and frequency domain
figure(1)
subplot(1,2,1) % time domain
plot(t,mex,'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('time (t)'); ylabel('\psi(t)')
axis([-6 6 -0.5 1.1])

subplot(1,2,2) % frequency domain
plot(ks,abs(fftshift(mext))/max(abs(mext)),'k','Linewidth',2)
set(gca,'Fontsize',16); xlabel('frequency (k)'); ylabel('FFT(\psi)')
axis([-10 10 0 1.05])




%% Create a function to return scaled and shifted Haar wavelets
function H = Haar(t,a,b)
    % INPUTS: 
    % t is the independent variable, input as a vector
    % a is the scale parameter
    % b is the translation parameter
    
    % OUTPUT:
    % H is the scaled and shifted Haar wavelet
    
    H = zeros(length(t),1);
    for j = 1:length(t)
        if t(j) >= b && t(j) < b + 0.5*a
            H(j) = 1/sqrt(abs(a));
        elseif t(j) >= b + 0.5*a && t(j) < b + a
            H(j) = -1/sqrt(abs(a));
        end
    end

end
