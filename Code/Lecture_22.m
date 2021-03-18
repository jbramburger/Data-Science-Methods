%% Lecture 22 Code

% Clean workspace
clear all; close all; clc

%% Create the 'A' touch-tone signal

N = 5000;
t = linspace(0,1/8,N);
f = sin(1394*pi*t) + sin(3266*pi*t);

figure(1)
plot(t,f,'Linewidth',2)
xlabel('t')
ylabel('f(t)')
set(gca,'Fontsize',16)
axis tight

%% Zoom in on signal

figure(2)
plot(t,f,'Linewidth',2)
xlabel('t')
ylabel('f(t)')
xlim([0 0.01])
set(gca,'Fontsize',16)

%% DCT of signal

ft = dct(f);

figure(3)
plot(ft,'r','Linewidth',2)
ylabel('DCT(f)')
set(gca,'Fontsize',16)

%% Zoom on peaks

figure(4)
plot(ft,'r','Linewidth',2)
ylabel('DCT(f)')
xlim([0 500])
set(gca,'Fontsize',16)

%% Sample 10% of the signal

m = 500; % number of points
perm = randperm(5000); % random permutation of the numbers from 1 to 5000
ind = perm(1:m); % 500 randomly selected numbers (indices of subsampled points)
tr = t(ind); % times of subsampled points
fr = f(ind); % subsampled values

figure(5)
plot(t,f,'b',tr,fr,'or','Linewidth',2)
xlim([0 0.01])
set(gca,'Fontsize',16)

%% Creating Psi and A matrices

Psi = idct(eye(N,N)); % Create matrix Psi that does inverse DCT
A = Psi(ind,:); % The A matrix is subsampled rows of Psi 

%% Solve for coefficients

b = fr'; % b is a column vector of subsampled values

x1 = A\b;
x2 = pinv(A)*b;

cvx_begin quiet
    variable x3(N);
    minimize( norm(x3,1) );
    subject to
        A*x3 == b;
cvx_end

%% Plot results

figure(6)
subplot(4,1,1)
plot(ft,'r','Linewidth',2)
xlim([0 500])
ylabel('DCT(f)')
title('exact')
set(gca,'Fontsize',16)
subplot(4,1,2)
plot(x1,'r','Linewidth',2)
xlim([0 500])
ylabel('x1')
title('backslash')
set(gca,'Fontsize',16)
subplot(4,1,3)
plot(x2,'r','Linewidth',2)
xlim([0 500])
ylabel('x2')
title('2-norm')
set(gca,'Fontsize',16)
subplot(4,1,4)
plot(x3,'r','Linewidth',2)
xlim([0 500])
ylabel('x3')
title('1-norm')
set(gca,'Fontsize',16)

%% Signal reconstruction with idct

sig1 = idct(x1);
sig2 = idct(x2);
sig3 = idct(x3);

figure(7)
subplot(4,1,1)
plot(t,f,'b','Linewidth',2)
ylabel('f')
xlabel('t')
title('exact')
set(gca,'Fontsize',16)
subplot(4,1,2)
plot(t,sig1,'b','Linewidth',2)
ylabel('f')
xlabel('t')
title('backslash')
set(gca,'Fontsize',16)
subplot(4,1,3)
plot(t,sig2,'b','Linewidth',2)
ylabel('f')
xlabel('t')
title('2-norm')
set(gca,'Fontsize',16)
subplot(4,1,4)
plot(t,sig3,'b','Linewidth',2)
ylabel('f')
xlabel('t')
title('1-norm')
set(gca,'Fontsize',16)

%% Zoom in on signals

figure(8)
subplot(4,1,1)
plot(t,f,'b',tr,fr,'or','Linewidth',2)
xlim([0 0.01])
ylabel('f')
xlabel('t')
title('exact')
set(gca,'Fontsize',16)
subplot(4,1,2)
plot(t,sig1,'b',tr,fr,'or','Linewidth',2)
xlim([0 0.01])
ylabel('f')
xlabel('t')
title('backslash')
set(gca,'Fontsize',16)
subplot(4,1,3)
plot(t,sig2,'b',tr,fr,'or','Linewidth',2)
xlim([0 0.01])
ylabel('f')
xlabel('t')
title('2-norm')
set(gca,'Fontsize',16)
subplot(4,1,4)
plot(t,sig3,'b',tr,fr,'or','Linewidth',2)
xlim([0 0.01])
ylabel('f')
xlabel('t')
title('1-norm')
set(gca,'Fontsize',16)