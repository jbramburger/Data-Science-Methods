%% Lecture 9 Code

% Clean workspace
clear all; close all; clc

%% Create the noisy dog image

I = imread('sherlock.jpg');
I = rgb2gray(I);
I = im2double(I);
In = imnoise(I,'gaussian',0,0.1);

figure(1)
imshow(I)

figure(2)
imshow(In)

%% Apply FFT to the noisy image

Int = fft2(In);
figure(3)
pcolor(fftshift(abs(Int)))
shading interp
axis([470 490 310 330])
colorbar

%% Plot FFT on a log scale

figure(3)
pcolor(log(fftshift(abs(Int))))
shading interp
axis([430 530 270 370])
colorbar

%% Creating the filter function

[Ny,Nx] = size(In);
kx = [0:Nx/2-1 -Nx/2:-1];
ky = [0:Ny/2-1 -Ny/2:-1];
[Kx,Ky] = meshgrid(kx,ky);

% Filter function
a = 1e-2;
filter = exp(-a*(Kx.^2 + Ky.^2));

%% Filtering the noisy image

Intf = Int.*filter;
Inf = ifft2(Intf);

figure(4)
subplot(2,2,1)
imshow(In)
subplot(2,2,2)
pcolor(fftshift(Kx),fftshift(Ky),log(fftshift(abs(Int))))
shading interp
subplot(2,2,3)
pcolor(fftshift(Kx),fftshift(Ky),log(fftshift(abs(filter))))
shading interp
subplot(2,2,4)
imshow(Inf)

%% Wider filter

a = 1e-3;
filter = exp(-a*(Kx.^2 + Ky.^2));

Intf = Int.*filter;
Inf = ifft2(Intf);

figure(5)
subplot(2,2,1)
imshow(In)
subplot(2,2,2)
pcolor(fftshift(Kx),fftshift(Ky),log(fftshift(abs(Int))))
shading interp
subplot(2,2,3)
pcolor(fftshift(Kx),fftshift(Ky),log(fftshift(abs(filter))))
shading interp
subplot(2,2,4)
imshow(Inf)


%% Comparing filter widths

% Different filtering widths
a = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7];

% Plot to compare filter widths
figure(6)
for j = 1:length(a)
   filter = exp(-a(j)*(Kx.^2 + Ky.^2));
   
   Intf = Int.*filter;
   Inf = ifft2(Intf);
   
   subplot(3,2,j)
   imshow(Inf)
   title(['a = ',num2str(a(j))])
end

%% Filtering with the Shannon filter

wx = 100;
wy = 100;
filter = zeros(size(Int));
filter([1:wy+1 Ny-wy+1:Ny],[1:wx+1 Nx-wx+1:Nx]) = ones(2*wy+1,2*wx+1);

Intf = Int.*filter;
Inf = ifft2(Intf);

figure(7)
subplot(2,2,1)
imshow(In)
subplot(2,2,2)
pcolor(fftshift(Kx),fftshift(Ky),log(fftshift(abs(Int))))
shading interp
subplot(2,2,3)
pcolor(fftshift(Kx),fftshift(Ky),log(fftshift(abs(filter))))
shading interp
subplot(2,2,4)
imshow(Inf)

%% Comparing different Shannon filter widths

% Different filtering widths
w = [25 50 100 150 200 250];

% Plot to compare filter widths
figure(8)
for j = 1:length(w)
   wx = w(j);
   wy = w(j);
   filter = zeros(size(Int));
   filter([1:wy+1 Ny-wy+1:Ny],[1:wx+1 Nx-wx+1:Nx]) = ones(2*wy+1,2*wx+1);

   
   Intf = Int.*filter;
   Inf = ifft2(Intf);
   
   subplot(3,2,j)
   imshow(Inf)
   title(['w = ',num2str(w(j))])
end

%% Mystery image to denoise

I = imread('mystery.jpg'); 
I = rgb2gray(I);
I = im2double(I);
In = imnoise(I,'gaussian',0,3);
figure(9)
imshow(In)

%% Shannon filtering of the mystery image

Int = fft2(In);
[Nx,Ny] = size(In);

% Different filtering widths
w = [25 50 100 150 200 250];

% Plot to compare filter widths
figure(10)
for j = 1:length(w)
   wx = w(j);
   wy = w(j);
   filter = zeros(size(Int));
   filter([1:wx+1 Nx-wx+1:Nx],[1:wy+1 Ny-wy+1:Ny]) = ones(2*wx+1,2*wy+1);

   Intf = Int.*filter;
   Inf = ifft2(Intf);
   
   subplot(3,2,j)
   imshow(Inf)
   title(['w = ',num2str(w(j))])
end
