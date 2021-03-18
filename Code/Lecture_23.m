%% Lecture 23 Code

% Clean workspace
clear all; close all; clc

%% Read in and print photo

A = imread('coffee.jpeg');
Abw2 = rgb2gray(A);
Abw = im2double(Abw2);
[nx, ny] = size(Abw);
figure(1), imshow(Abw)

%% FFT of image

At = fftshift(fft2(Abw));

figure(2)
subplot(2,1,1)
plot(reshape(abs(At),nx*ny,1),'k','Linewidth',2)
ylabel('|ut|')
set(gca,'Xlim',[2.3*10^5 2.5*10^5],'Fontsize',16)

% semi-log plot
subplot(2,1,2)
semilogy(reshape(abs(At),nx*ny,1),'k','Linewidth',2)
ylabel('log(|ut|)')
set(gca,'Xlim',[2.3*10^5 2.5*10^5],'Fontsize',16)

%% Threshold Fourier coefficients (takes a few minutes)

% Open new figure
figure(3)

% Plot original image
subplot(2,2,1)
imshow(Abw)

% Plot thresholded images
count_pic = 2;
for thresh = [1e2 5e2 1e3] 
    At2 = reshape(At,nx*ny,1);
    count = 0;
    for j = 1:length(At2)
       if abs(At2(j)) < thresh
          At2(j) = 0;
          count = count + 1;
       end
    end
    
    percent = 100 - (count/length(At2))*100
    
    Atlow = fftshift(reshape(At2,nx,ny)); 
    Alow = (ifft2(Atlow));
    subplot(2,2,count_pic)
    imshow(Alow)
    count_pic = count_pic + 1;
end

%% Resize image to speed up computations

B = imresize(A,[75,100]);

Abw2 = rgb2gray(B);
Abw = im2double(Abw2);
[ny, nx] = size(Abw);

%% Sample the image

k = 2500; % number of sparse samples
Atest2 = zeros(ny,nx);
perm = randperm(nx*ny);
ind = perm(1:k);
for j = 1:k
   Atest2(ind(j)) = -1; 
end
imshow(Atest2), caxis([-1 0])

%% Initialize underdetermined system (takes a minute)

Atest = zeros(ny,nx);
for j = 1:k
   Atest(ind(j)) = 1;
   Adel = reshape(idct2(Atest),nx*ny,1);
   Adelta(j,:) = Adel;
   Atest(ind(j)) = 0;
end

b = Abw(ind)';

%% L^1 Optimization (takes a while to run)

n = nx*ny;
cvx_begin quiet
    variable y(n);
    minimize( norm(y,1) );
    subject to
        Adelta*y == b;
cvx_end

Alow = dct2(reshape(y,ny,nx));

figure(4)
subplot(1,2,1)
imshow(Abw)
subplot(1,2,2)
imshow(Alow)





