%% Lecture 12 Code

% Clean workspace
clear all; close all; clc


%% Example 1 function

x = linspace(0,1,25);
t = linspace(0,2,50);
[T, X] = meshgrid(t,x);
f = exp(-abs((X - 0.5).*(T - 1))) + sin(X.*T);
figure(1)
surf(X,T,f)
shading interp
set(gca,'fontsize',16)

%% Example 1 rank-1 approximation

[U,S,V] = svd(f,'econ');
f_rank1 = U(:,1)*S(1,1)*V(:,1)';

figure(2)
subplot(1,2,1)
surf(X,T,f)
shading interp
subplot(1,2,2)
surf(X,T,f_rank1)
shading interp

%% Example 1 comparing low-rank approximations

figure(3)
subplot(2,2,1)
surf(X,T,f), shading interp
for j = 1:3
   ff = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
   subplot(2,2,j+1)
   surf(X,T,ff), shading interp
   set(gca,'Zlim',[0.5 2])
end

%% Calculating the energy of the truncations

sig = diag(S);
energy1 = sig(1)^2/sum(sig.^2)
energy2 = sum(sig(1:2).^2)/sum(sig.^2)
energy3 = sum(sig(1:3).^2)/sum(sig.^2)

%% Energies captured by each mode

figure(4)
subplot(2,2,1)
plot(sig,'ko','Linewidth',2)
axis([0 25 0 50])
ylabel('\sigma')
set(gca,'Fontsize',16,'Xtick',0:5:25)
subplot(2,2,2)
semilogy(sig,'ko','Linewidth',2)
axis([0 25 10^(-18) 10^5])
ylabel('\sigma (log scale)')
set(gca,'Fontsize',16,'Xtick',0:5:25,'Ytick',logspace(-15,5,5))
subplot(2,2,3)
plot(sig.^2/sum(sig.^2),'ko','Linewidth',2)
axis([0 25 0 1])
ylabel('Energy')
set(gca,'Fontsize',16,'Xtick',0:5:25)
subplot(2,2,4)
semilogy(sig.^2/sum(sig.^2),'ko','Linewidth',2)
axis([0 25 10^(-18) 10^5])
ylabel('Energy (log scale)')
set(gca,'Fontsize',16,'Xtick',0:5:25,'Ytick',logspace(-15,0,4))

%% Cumulative energy in the first N modes

figure(5)
subplot(1,2,1)
plot(cumsum(sig.^2)/sum(sig.^2),'ko','Linewidth',2)
axis([0 25 10^-(18) 1])
ylabel('Cumulative Energy')
set(gca,'Fontsize',16,'Xtick',0:5:25)
subplot(1,2,2)
semilogy(cumsum(sig.^2)/sum(sig.^2),'ko','Linewidth',2)
axis([0 25 10^-(18) 1])
ylabel('Cumulative Energy (log scale)')
set(gca,'Fontsize',16,'Xtick',0:5:25,'Ytick',logspace(-15,0,4))

%% Plotting the POD modes

figure(6)
subplot(3,2,1)
surf(X,T,f), shading interp
xlabel('x')
ylabel('t')
title('f','Fontsize',16)
set(gca,'Fontsize',16)
subplot(3,2,2)
surf(X,T,f_rank1), shading interp
xlabel('x')
ylabel('t')
title('rank 1','Fontsize',16)
set(gca,'Fontsize',16)
subplot(3,1,2)
plot(x,U(:,1),'b',x,U(:,2),'--r',x,U(:,3),':k','Linewidth',2)
xlabel('x')
set(gca,'Fontsize',16)
legend('mode 1','mode 2','mode 3','Location','best')
subplot(3,1,3)
plot(t,V(:,1),'b',t,V(:,2),'--r',t,V(:,3),':k','Linewidth',2)
xlabel('t')
set(gca,'Fontsize',16)


%% Example 2 waterfall plot

clear all; close all; clc

x = linspace(-10,10,100);
t = linspace(0,10,30);
[X,T] = meshgrid(x,t);
f = sech(X).*(1-0.5*cos(2*T))+(sech(X).*tanh(X)).*(1-0.5*sin(2*T));
figure(7)
waterfall(X,T,f)
xlabel('x')
ylabel('t')
set(gca,'Fontsize',16,'Zlim',[-1,2])

%% Plotting the low-rank approximations for Example 2

[U,S,V] = svd(f','econ');
figure(8)
subplot(2,2,1)
waterfall(X,T,f)
xlabel('x')
ylabel('t')
set(gca,'Fontsize',16,'Zlim',[-1,2])
for j=1:3
    ff = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
    subplot(2,2,j+1)
    waterfall(X,T,ff')
    set(gca,'Fontsize',16,'Zlim',[-1,2])
end

%% Plotting the POD modes for Example 2

sig = diag(S);
figure(9)
subplot(3,2,1)
plot(sig,'ko','Linewidth',2)
axis([0 25 0 50])
ylabel('\sigma')
set(gca,'Fontsize',16,'Xtick',0:5:25)
subplot(3,2,2)
semilogy(sig,'ko','Linewidth',2)
axis([0 25 10^-(18) 10^5])
ylabel('\sigma (log scale)')
set(gca,'Fontsize',16,'Xtick',0:5:25,'Ytick',logspace(-15,5,5))
subplot(3,1,2)
plot(x,U(:,1),'b',x,U(:,2),'--r',x,U(:,3),':k','Linewidth',2)
set(gca,'Fontsize',16)
legend('mode 1','mode 2','mode 3','Location','northwest')
subplot(3,1,3)
plot(t,V(:,1),'b',t,V(:,2),'--r',t,V(:,3),':k','Linewidth',2)
set(gca,'Fontsize',16)
legend('v_1','v_2','v_3','Location','northwest')