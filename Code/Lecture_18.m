%% Lecture 18 Code

% Clean workspace
clear all; close all; clc

%% Loading cat and dog data

load('catData.mat')
load('dogData.mat')

%% Plotting first 9 dog images

figure(1)
for k = 1:9
   subplot(3,3,k)
   dog1 = reshape(dog(:,k),64,64);
   imshow(dog1)
end

%% Discrete wavelet transform of dog image

X = im2double(reshape(dog(:,6),64,64));
[cA, cH, cV, cD] = dwt2(X,'haar');

%% Plotting dwt results

figure(2)
subplot(2,2,1)
imshow(cA)
subplot(2,2,2)
imshow(cH)
subplot(2,2,3)
imshow(cV)
subplot(2,2,4)
imshow(cD)

%% Rescaling the results of the dwt

cod_cH1 = rescale(abs(cH));
cod_cV1 = rescale(abs(cV));
cod_edge = cod_cH1+cod_cV1;
figure(3)
subplot(2,2,1)
imshow(cod_cH1)
subplot(2,2,2)
imshow(cod_cV1)
subplot(2,2,3)
imshow(cod_edge)
subplot(2,2,4)
imshow(X)
