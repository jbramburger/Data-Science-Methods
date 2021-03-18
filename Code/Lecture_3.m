%% Lecture 3 Code

% Clean workspace
clear all; close all; clc


%% Creating a noisy signal

L = 30; % timeslot [-L,L]
n = 512; % number of Fourier modes

t2 = linspace(-L,L,n+1);
t = t2(1:n);
k = (2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; %define this so we don't need to use fftshift(k) every time we plot
ks = fftshift(k);
u = sech(t);
ut = fft(u);

% Add noise
noise = 10;
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));

% Plot noisy signal
figure(1)
plot(ks,fftshift(abs(utn))/max(abs(utn)),'r','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('frequency (k)')
ylabel('|ut|')

%% Average over 2 realizations of the signal

ave = zeros(1,n);
for j = 1:2
   utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
   ave = ave + utn;
end
ave = abs(fftshift(ave))/2;

% Plot averaged signal
figure(2)
plot(ks,ave/max((ave)),'r','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('frequency (k)')
ylabel('|ut|')

%% Average over 5 realizations of the signal

ave = zeros(1,n);
for j = 1:5
   utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
   ave = ave + utn;
end
ave = abs(fftshift(ave))/5;

% Plot averaged signal
figure(3)
plot(ks,ave/max((ave)),'r','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('frequency (k)')
ylabel('|ut|')

%% Comparing over different realizations

labels = ['(a)';'(b)';'(c)';'(d)'];
realize = [1 2 5 200];

figure(4)
for jj = 1:length(realize)
   u = sech(t);
   ave = zeros(1,n);
   ut = fft(u);
   
   % Averaging over realizations
   for j = 1:realize(jj)
      utn(j,:) = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
      ave = ave+utn(j,:);
   end
   ave = abs(fftshift(ave))/realize(jj);
   
   % Plot averaged signals in one plot
   subplot(length(realize),1,jj)
   plot(ks,ave/max(ave),'r','Linewidth',2)
   set(gca,'Fontsize',16)
   axis([-20 20 0 1])
   text(-18,0.7,labels(jj,:),'Fontsize',16)
   ylabel('|fft(u)|','Fontsize',16)
end

hold on
plot(ks,abs(fftshift(ut))/max(abs(ut)),'k:','Linewidth',2)
xlabel('frequency (k)')

%% Shifted realizations

slice = 0:0.5:10;
[T,S] = meshgrid(t,slice);
[K,S] = meshgrid(k,slice);

U = sech(T - 10*sin(S)).*exp(1i*0*T);
figure(5)
subplot(2,1,1)
waterfall(T,2*S+1,U), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-30 30],'Zlim',[0 2])
xlabel('time (t)'), ylabel('realizations'), zlabel('|u|')

for j = 1:length(slice)
   Ut(j,:) = fft(U(j,:));
   Kp(j,:) = fftshift(K(j,:));
   Utp(j,:) = fftshift(Ut(j,:));
   Utn(j,:) = Ut(j,:) + + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
   Utnp(j,:) = fftshift(Utn(j,:))/max(abs(Utn(j,:)));
   Un(j,:) = ifft(Utn(j,:));
end

subplot(2,1,2)
waterfall(Kp,2*S+1,abs(Utp)/max(abs(Utp(1,:)))), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-28 28])
xlabel('frequency (k)'), ylabel('realizations'), zlabel('|fft(u)|')

%% Noisy realizations

figure(6)
subplot(2,1,1)
waterfall(T,2*S+1,abs(Un)), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-30 30],'Zlim',[0 2])
xlabel('time (t)'), ylabel('realizations'), zlabel('|u|')

subplot(2,1,2)
waterfall(Kp,2*S+1,abs(Utnp)/max(abs(Utnp(1,:)))), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-28 28])
xlabel('time (t)'), ylabel('realizations'), zlabel('|fft(u)|')

%% Average over time and frequency domains to compare

Uave = zeros(1,n);
Utave = zeros(1,n);
for j = 1:length(slice)
   Uave = Uave + Un(j,:);
   Utave = Utave + Utn(j,:);
end
Uave = Uave/length(slice);
Utave = fftshift(Utave)/length(slice);

figure(7) 
subplot(2,1,1)
plot(t,abs(Uave),'k')
set(gca,'Fontsize',16)
xlabel('time (t)'), ylabel('|u|')

subplot(2,1,2)
plot(t,abs(Utave)/max(abs(Utave)),'k')
hold on
plot(ks,abs(fftshift(Ut(1,:))/max(abs(Ut(1,:)))),'k','Linewidth',2)
axis([-20 20 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('|fft(u)|')
