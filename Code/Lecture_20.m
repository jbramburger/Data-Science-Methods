%% Lecture 20 Code

% Clean workspace
clear all; close all; clc

%% Train model

load('catData.mat')
load('dogData.mat')

dog_wave = dc_wavelet(dog);
cat_wave = dc_wavelet(cat);

feature = 20;
[U,S,V,threshold,w,sortdog,sortcat] = dc_trainer(dog_wave,cat_wave,feature);

%% Load and view from the test set

load('PatternRecAns')

figure(1)
for k = 1:9
   subplot(3,3,k)
   test = reshape(TestSet(:,k+5),64,64);
   imshow(test)
end

%% Classify test data

TestNum = size(TestSet,2);
Test_wave = dc_wavelet(TestSet); % wavelet transform
TestMat = U'*Test_wave; % PCA projection
pval = w'*TestMat;

%% Check pval against threshold

% Cat = 1, dog = 0
ResVec = (pval > threshold)

%% Checking performance

% 0s are correct and 1s are incorrect
err = abs(ResVec - hiddenlabels)
errNum = sum(err);
sucRate = 1 - errNum/TestNum;

%% Checking errors

k = 1;
TestNum = length(pval);
figure(2)
for j = 1:TestNum
   if ResVec(j) ~= hiddenlabels(j)
      S = reshape(TestSet(:,j),64,64);
      subplot(1,2,k)
      imshow(S)
      k = k+1;
   end
end

%% Classifying your pets

% Clear pval
pval = [];

close all

% Read in image
%I = imread('Bailey.jpg'); % Change the file name for your photo
%I = imread('Lola.jpg');
%I = imread('Dingbat.jpg');
%I = imread('Carter_Peyton_Oakey.png');
%I = imread('Zach_Zlepper_Duke.jpg');
%I = imread('Eakin_Shen_Shio.jpg');
%I = imread('Eakin_Shen_Shio_2.jpg');
%I = imread('Haley_Riggs_AJ.jpg');
%I = imread('Haley_Riggs_Latte.jpg');
%I = imread('Christian_Valoria_Luna.png');
%I = imread('Aidan_Hunt_Io.jpeg');
%I = imread('Alex_Troy_Mallen_Prince.jpeg');
%I = imread('Haoran_Sun_Rich.jpeg');
%I = imread('John_Curtis_Cody.jpg');
%I = imread('Sophia_Jannetty_Chicky.jpg');
imshow(I)

% Convert to grayscale and resize
I = rgb2gray(I);
I = im2double(I);
I = imresize(I,[64,64]);
I = reshape(I,64*64,1);

%  Classify the image
I_wave = dc_wavelet(I); % wavelet transform
IMat = U'*I_wave; % PCA projection
pval = w'*IMat;

if pval > threshold
    disp('Cat')
else
    disp('Dog')
end
