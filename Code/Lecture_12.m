%% Lecture 12 Code

% Clean workspace
clear all; close all; clc


%% Load and plot the weight and height data

load('weightheight.mat')

plot(X(1,:),X(2,:),'k.','MarkerSize',10)
hold on
axis equal
xlabel('weight')
ylabel('height')
set(gca,'Fontsize',16)
















%% Calculate the SVD of the weight-height matrix

[U,S,V] = svd(X,'econ');



















%% Rank-1 Approximation of weight-height matrix

X_rank1 = S(1,1)*U(:,1)*V(:,1)';
plot(X_rank1(1,:),X_rank1(2,:),'r.','MarkerSize',10)




















%% Projecting onto the rank-1 approximation

plot([X(1,10) X_rank1(1,10)],[X(2,10) X_rank1(2,10)],'g','Linewidth',2)
























%% Plotting the line for the rank-1 approximation

n = 200;
y1 = S(1,1)/sqrt(n-1)*U(:,1);
c = compass(y1(1),y1(2)); % creates vector from the origin pointed at y1
set(c,'Linewidth',4)



























%% Plotting the line in the direction of u_2

y2 = S(2,2)/sqrt(n-1)*U(:,2);
c = compass(y2(1),y2(2)); 
set(c,'Linewidth',4)





















%% Change of basis using U

X_proj = U'*X;

figure(2)
plot(X_proj(1,:),X_proj(2,:),'k.','MarkerSize',10)
axis equal
hold on
y1_proj = U'*y1;
y2_proj = U'*y2;
c = compass(y1_proj(1),y1_proj(2));
set(c,'Linewidth',4);
c = compass(y2_proj(1),y2_proj(2));
set(c,'Linewidth',4);





















%% Low-rank approximations of 3D data

% Clean workspace
clear all; close all; clc

rng(5); % make random numbers be the same every time
n = 200; % number of data points
X = [3 6 4; 6 3 0; 4 0 1]*randn(3,n); % create n data points in 3D

[U,S,V] = svd(X,'econ');

figure(3)
plot3(X(1,:),X(2,:),X(3,:),'k.','Markersize',10)
axis vis3d
hold on
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'Fontsize',16)

% Compute the rank-1 approximation of the data
X_rank1 = U(:,1)*S(1,1)*V(:,1)';
plot3(X_rank1(1,:),X_rank1(2,:),X_rank1(3,:),'r.','Markersize',10)
 
% Plot the vector sqrt(sigma_1)*u_1/(n-1)
vec = S(1,1)/sqrt(n-1) * U(:,1);
quiver3(0,0,0, vec(1),vec(2),vec(3),0,'b','Linewidth',5,'maxheadsize',1)
























%% Vector point in u_2 direction
vec = S(2,2)/sqrt(n-1) * U(:,2);
quiver3(0,0,0, vec(1),vec(2),vec(3),0,'b','Linewidth',5,'maxheadsize',1)


















%% Rank-2 approximation

X_rank2 = U(:,1)*S(1,1)*V(:,1)' + U(:,2)*S(2,2)*V(:,2)';
plot3(X_rank2(1,:),X_rank2(2,:),X_rank2(3,:),'g.','Markersize',10)























%% Another 3D dataset

% Clean workspace
clear all; close all; clc

t = -5:0.05:5;
x = 5*t;
y = 4*t;
z = -6*t;
X = [x; y; z];

figure(4)
plot3(X(1,:),X(2,:),X(3,:),'k.','Markersize',10)
axis vis3d
hold on
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'Fontsize',16)


















%% Calculate the SVD

[U,S,V] = svd(X,'econ');






















%% Rank-1 Approximation

X_rank1 = U(:,1)*S(1,1)*V(:,1)';
plot3(X_rank1(1,:),X_rank1(2,:),X_rank1(3,:),'r.','Markersize',10)


















