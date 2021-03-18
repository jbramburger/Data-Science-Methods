%% Lecture 8 Code

% Clean workspace
clear all; close all; clc












%% Load and show an image
I = imread('sherlock.jpg');
imshow(I)




























%% Convert to grayscale
Ig = rgb2gray(I);
imshow(Ig)




















%% Add noise
noiseg = randn(640,960);
Ign = Ig + noiseg;
















%% Convert the image to double precision and add noise

I = im2double(I);
Ig = im2double(Ig);
Ign = Ig + noiseg;
imshow(Ign)

noise = randn(640,960,3);
In = I + noise;
imshow(In)





















%% MATLAB routine for adding noise

Ign = imnoise(Ig,'gaussian');
In = imnoise(I,'gaussian');




















%% Plotting all the images

subplot(2,2,1)
imshow(I)
subplot(2,2,2)
imshow(In)
subplot(2,2,3)
imshow(Ig)
subplot(2,2,4)
imshow(Ign)




















%% Changing the mean and variance of the noise

Ign = imnoise(Ig,'gaussian',0,0.1);
In = imnoise(I,'gaussian',0.1,0.01);

subplot(2,2,1)
imshow(I)
subplot(2,2,2)
imshow(In)
subplot(2,2,3)
imshow(Ig)
subplot(2,2,4)
imshow(Ign)



















%% Salt and pepper noise

Ign = imnoise(Ig,'salt & pepper',0.5);
In = imnoise(I,'salt & pepper',0.5);

subplot(2,2,1)
imshow(I)
subplot(2,2,2)
imshow(In)
subplot(2,2,3)
imshow(Ig)
subplot(2,2,4)
imshow(Ign)
