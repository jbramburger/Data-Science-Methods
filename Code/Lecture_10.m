%% Lecture 10 Code

% Clean workspace
clear all; close all; clc


%% Create the clean and noisy signal

L = 30; n =128;
x2 = linspace(-L,L,n+1);
x = x2(1:n);

u0 = sech(x); 
ut = fft(u0);
noise = 2;
utn = ut + noise*(randn(1,n) + 1i*randn(1,n));
un = real(ifft(utn));

figure(1)
pos = get(gcf,'Position');
set(gcf,'Position',[0 100 3*pos(3) pos(4)])
subplot(1,3,1)
plot(x,u0,'Linewidth',2)
axis([-L L -0.5 1.2])
xlabel('x','Fontsize',16)
title('Clean Signal','Fontsize',16)
subplot(1,3,2)
plot(x,un,'Linewidth',2)
axis([-L L -0.5 1.2])
xlabel('x','Fontsize',16)
title('Noisy Signal','Fontsize',16)

%% Heat equation solver

D = 1;
h = x(2) - x(1);

% Create A matrix
e = ones(n,1);
A = spdiags([e -2*e e], [-1 0 1],n,n);
A(1,n) = 1; % periodic BCs
A(n,1) = 1;

% ode45 for time-stepping
rhs = @ (t,u) (D/h^2)*A*u;
[t,u] = ode45(rhs,linspace(0,0.5,101),un);

% Plot a video of the evolving heat equation solution
figure(1)
for j = 1:length(t)
   subplot(1,3,3)
   plot(x,u(j,:),'Linewidth',2)
   axis([-L L -0.5 1.2])
   xlabel('x','Fontsize',16)
   title(['Denoised Signal, t = ',num2str(t(j))],'Fontsize',16)
   pause(0.1)
end

%% Creating a noisy image 

I = imread('cameraman.tif');
I = im2double(I);
In = imnoise(I,'gaussian',0,0.1);

figure(2)
pos = get(gcf,'Position');
set(gcf,'Position', [0 100 3*pos(3) pos(4)])
subplot(1,3,1)
imshow(I)
title('Clean Image','Fontsize',16)
subplot(1,3,2)
imshow(In)
title('Noisy Image','Fontsize',16)

%% Denoising with the 2D heat equation

D = 1;
[ny,nx] = size(I);
x = linspace(0,1,nx);
y = linspace(0,1,ny);
dx = x(2) - x(1);
dy = y(2) - y(1);
[X,Y] = meshgrid(x,y);

% Build matrix A using sparse matrices
ex = ones(nx,1);
ey = ones(ny,1);
Dx = (spdiags([ex -2*ex ex],[-1 0 1],nx,nx))/dx^2;
Ix = speye(nx);
Dy = (spdiags([ey -2*ey ey],[-1 0 1],ny,ny))/dy^2;
Iy = speye(ny);
L = kron(Iy,Dx) + kron(Dy,Ix);

% Simulate ODE
rhs = @ (t,u) D*L*u;
[t,u] = ode45(rhs,linspace(0,2.5e-5,101),reshape(In,nx*ny,1));

% Plot denoised image
figure(2)
for j = 1:length(t)
   subplot(1,3,3)
   U = reshape(u(j,:).',ny,nx);
   imshow(U)
   drawnow
   title(['Denoised Image, t = ',num2str(t(j))],'Fontsize',16)
   pause(0.1)
end




