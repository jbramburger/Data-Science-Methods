%% Lecture 19 Code

% Clean workspace
clear all; close all; clc

%% Load in data and apply DWT

load('catData.mat')
load('dogData.mat')

dog_wave = dc_wavelet(dog);
cat_wave = dc_wavelet(cat);

%% Apply SVD to dog and cat data

[U,S,V] = svd([dog_wave cat_wave],'econ');

%% Plot first four principal components

for k = 1:4
   subplot(2,2,k)
   ut1 = reshape(U(:,k),32,32);
   ut2 = rescale(ut1);
   imshow(ut2)
end

%% Plot singular values

figure(2)
subplot(2,1,1)
plot(diag(S),'ko','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 80])
subplot(2,1,2)
semilogy(diag(S),'ko','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 80])

%% Plot right singular vectors

figure(3) 
for k = 1:3
   subplot(3,2,2*k-1)
   plot(1:40,V(1:40,k),'ko-')
   subplot(3,2,2*k)
   plot(1:40,V(81:120,k),'ko-')
end
subplot(3,2,1), set(gca,'Ylim',[-.15 0],'Fontsize',12), title('dogs')
subplot(3,2,2), set(gca,'Ylim',[-.15 0],'Fontsize',12), title('cats')
subplot(3,2,3), set(gca,'Ylim',[-.2 .2],'Fontsize',12)
subplot(3,2,4), set(gca,'Ylim',[-.2 .2],'Fontsize',12)
subplot(3,2,5), set(gca,'Ylim',[-.2 .2],'Fontsize',12)
subplot(3,2,6), set(gca,'Ylim',[-.2 .2],'Fontsize',12)

%% Project onto PCA modes

feature = 20;

nd = size(dog_wave,2);
nc = size(cat_wave,2);
animals = S*V'; % projection onto principal components: X = USV' --> U'X = SV'
dogs = animals(1:feature,1:nd);
cats = animals(1:feature,nd+1:nd+nc);

%% Calculate scatter matrices

md = mean(dogs,2);
mc = mean(cats,2);

Sw = 0; % within class variances
for k = 1:nd
    Sw = Sw + (dogs(:,k) - md)*(dogs(:,k) - md)';
end
for k = 1:nc
   Sw =  Sw + (cats(:,k) - mc)*(cats(:,k) - mc)';
end

Sb = (md-mc)*(md-mc)'; % between class

%% Find the best projection line

[V2, D] = eig(Sb,Sw); % linear disciminant analysis
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

%% Project onto w

vdog = w'*dogs;
vcat = w'*cats;

%% Make dogs below the threshold

if mean(vdog) > mean(vcat)
    w = -w;
    vdog = -vdog;
    vcat = -vcat;
end

%% Plot dog/cat projections (not for function)

figure(4)
plot(vdog,zeros(80),'ob','Linewidth',2)
hold on
plot(vcat,ones(80),'dr','Linewidth',2)
ylim([0 1.2])

%% Find the threshold value

sortdog = sort(vdog);
sortcat = sort(vcat);

t1 = length(sortdog);
t2 = 1;
while sortdog(t1) > sortcat(t2)
    t1 = t1 - 1;
    t2 = t2 + 1;
end
threshold = (sortdog(t1) + sortcat(t2))/2;

%% Plot histogram of results

figure(5)
subplot(1,2,1)
histogram(sortdog,30); hold on, plot([threshold threshold], [0 10],'r')
set(gca,'Xlim',[-3 4],'Ylim',[0 10],'Fontsize',14)
title('dog')
subplot(1,2,2)
histogram(sortcat,30); hold on, plot([threshold threshold], [0 10],'r')
set(gca,'Xlim',[-3 4],'Ylim',[0 10],'Fontsize',14)
title('cat')