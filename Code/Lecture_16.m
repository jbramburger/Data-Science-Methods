%% Lecture 16 Code

% Clean workspace
clear all; close all; clc


%% Load and view image files

S1 = imread('sicily1','jpeg');
S2 = imread('sicily2','jpeg');

subplot(1,2,1), imshow(S1);
subplot(1,2,2), imshow(S2);



%% Mixing the images 

A = [4/5 3/5; 1/2 2/3];
X1 = double(A(1,1)*S1 + A(1,2)*S2);
X2 = double(A(2,1)*S1 + A(2,2)*S2);

subplot(1,2,1), imshow(uint8(X1));
subplot(1,2,2), imshow(uint8(X2));


%% Create U^* matrix

% Reshape data
[m,n] = size(X1);
x1 = reshape(X1,m*n,1);
x2 = reshape(X2,m*n,1);

% Remove mean
x1 = x1 - mean(x1); x2 = x2 - mean(x2); 

% U^* matrix
theta0 = 0.5*atan(-2*sum(x1.*x2)/sum(x1.^2 - x2.^2));
Us = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];


%% Finding Sigma^{-1} matrix

sig1 = sum( (x1*cos(theta0) + x2*sin(theta0)).^2);
sig2 = sum( (x1*cos(theta0 - pi/2) + x2*sin(theta0 - pi/2)).^2);
Sigma = [1/sqrt(sig1) 0; 0 1/sqrt(sig2)];


%% Creating V matrix

% Create barred variables
X1bar = Sigma(1,1)*(Us(1,1)*X1 + Us(1,2)*X2);
X2bar = Sigma(2,2)*(Us(2,1)*X1 + Us(2,2)*X2);

x1bar = reshape(X1bar,m*n,1);
x2bar = reshape(X2bar,m*n,1);

% V matrix
phi0 = 0.25*atan( -sum( (2*(x1bar.^3).*x2bar - 2*x1bar.*(x2bar.^3)))... 
    /sum(3*(x1bar.^2).*(x2bar.^2) - 0.5*(x1bar.^4) - 0.5*(x2bar.^4)) );

V = [cos(phi0) sin(phi0); -sin(phi0) cos(phi0)];

% Recover independent components
S1bar = V(1,1)*X1bar + V(1,2)*X2bar;
S2bar = V(2,1)*X1bar + V(2,2)*X2bar;


%% Rescale and plot final images

min1 = min(min(min(S1bar)));
S1bar = S1bar - min1;
max1 =  max(max(max(S1bar)));
S1bar = S1bar*(255/max1);

min2 = min(min(min(S2bar)));
S2bar = S2bar - min2;
max2 =  max(max(max(S2bar)));
S2bar = S2bar*(255/max2);

% Plot mixed images and the independent components together

subplot(2,2,1), imshow(uint8(X1));
subplot(2,2,2), imshow(uint8(X2));
subplot(2,2,3), imshow(uint8(S1bar));
subplot(2,2,4), imshow(uint8(S2bar));

%% Mixing songs

% Clean workspace
clear all; close all; clc

load('laughter')
S1 = y;
load('handel')
S2 = y(1:length(S1));

A = [0.85 -0.5; 0.3 0.9];
X1 = double(A(1,1)*S1 + A(1,2)*S2);
X2 = double(A(2,1)*S1 + A(2,2)*S2);

% Play the mixed signals
% soundsc(X1)
% soundsc(X2)

% Remove mean
x1 = X1 - mean(X1); x2 = X2 - mean(X2); 

% U^* matrix 
theta0 = 0.5*atan(-2*sum(x1.*x2)/sum(x1.^2 - x2.^2));
Us = [cos(theta0) -sin(theta0); sin(theta0) cos(theta0)];

% Sigma values
sig1 = sum( (x1*cos(theta0) + x2*sin(theta0)).^2);
sig2 = sum( (x1*cos(theta0 - pi/2) + x2*sin(theta0 - pi/2)).^2);
Sigma = [1/sqrt(sig1) 0; 0 1/sqrt(sig2)];

% Create barred variables
x1bar = Sigma(1,1)*(Us(1,1)*X1 + Us(1,2)*X2);
x2bar = Sigma(2,2)*(Us(2,1)*X1 + Us(2,2)*X2);

% V matrix
phi0 = 0.25*atan( -sum( (2*(x1bar.^3).*x2bar - 2*x1bar.*(x2bar.^3)))... 
    /sum(3*(x1bar.^2).*(x2bar.^2) - 0.5*(x1bar.^4) - 0.5*(x2bar.^4)) );

V = [cos(phi0) sin(phi0); -sin(phi0) cos(phi0)];

% Recover independent components
S1bar = V(1,1)*x1bar + V(1,2)*x2bar;
S2bar = V(2,1)*x1bar + V(2,2)*x2bar;







