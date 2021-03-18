%% Lecture 2 Code

% Clean workspace
clear all; close all; clc


%% Fourier transform of sech(t)

L = 30;	% time slot to transform
n = 512;	% number of Fourier modes 2^9
t2 = linspace(-L,L,n+1);	% time discretization 
t = t2(1:n);		% only use the first n points (periodicity)
k = (2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 	%frequency components
u = sech(t);	%ideal signal in the time domain

% Plot the time domain
figure(1)
subplot(2,1,1)
plot(t,u, 'Linewidth', 2)
set(gca, 'Fontsize', 16)
xlabel('time (t)')
ylabel('|u|')

% Plot the frequency domain
ut = fft(u);
subplot(2,1,2)
plot(fftshift(k), fftshift(abs(ut)), 'r', 'Linewidth', 2)
set(gca, 'Fontsize', 16)
xlabel('frequency (k)')
ylabel('|ut|')

%% Adding noise

% Plot ideal signal
figure(2)
subplot(3,1,1)
plot(t,u,'Linewidth',2)
xlabel('time (t)')
ylabel('u')
set(gca,'Fontsize',16)

% Plot noisy signal (noise amplitude = 1)
noise = 1;
% Add white noise to the signal in frequency domain
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
un = ifftshift(utn); % Noisy signal in time domain
subplot(3,1,2)
plot(t,abs(un),'Linewidth',2)
xlabel('time (t)')
ylabel('u')
set(gca,'Fontsize',16)

% Plot noisy signal (noise amplitude = 5)
noise = 5;
% Add white noise to the signal in frequency domain
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
un = ifftshift(utn); % Noisy signal in time domain
subplot(3,1,3)
plot(t,abs(un),'Linewidth',2)
xlabel('time (t)')
ylabel('u')
set(gca,'Fontsize',16)

%% Even noisier signal

noise = 10;
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
un = ifft(utn);

% Plot the noisy signal 
figure(3)
subplot(2,1,1)
plot(t,abs(un),'Linewidth',2)
axis([-30 30 0 2])
xlabel('time (t)')
ylabel('|u|')
set(gca,'Fontsize',16)

% Plot the Fourier transform
subplot(2,1,2)
plot(fftshift(k), abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')
set(gca,'Fontsize',16)

%% Frequency filtering with a Gaussian

tau = 0.2;
k0 = 0;
filter = exp(-tau*(k - k0).^2); % Define the filter
unft = filter.*utn; % Apply the filter to the signal in frequency space
unf = ifft(unft);

% Plot the unfiltered signal in the frequency domain and the Gaussian filter
figure(4)
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the frequency domain
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the time domain and the ideal signal
subplot(3,1,3)
plot(t,u,'k','Linewidth',2)
hold on
plot(t,abs(unf),'Linewidth',2)
xlabel('time (t)')
ylabel('|u|')

%% No signal to filter

utn = noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
tau = 0.2;
k0 = 0;
filter = exp(-tau*(k - k0).^2); % Define the filter
unft = filter.*utn; % Apply the filter to the signal in frequency space
unf = ifft(unft);

% Plot the unfiltered signal in the frequency domain and the Gaussian
% filter
figure(5)
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the frequency domain
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the time domain and the ideal signal
subplot(3,1,3)
plot(t,u,'k','Linewidth',2)
hold on
plot(t,abs(unf),'Linewidth',2)
xlabel('time (t)')
ylabel('|u|')

%% Frequency filtering around the wrong frequency

tau = 0.2;
k0 = 15;
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
filter = exp(-tau*(k - k0).^2); % Define the filter
unft = filter.*utn; % Apply the filter to the signal in frequency space
unf = ifft(unft);

% Plot the unfiltered signal in the frequency domain and the Gaussian
% filter
figure(6)
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the frequency domain
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the time domain and the ideal signal
subplot(3,1,3)
plot(t,u,'k','Linewidth',2)
hold on
plot(t,abs(unf),'Linewidth',2)
xlabel('time (t)')
ylabel('|u|')