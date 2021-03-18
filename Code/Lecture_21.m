%% Lecture 21 Code

% Clean workspace
clear all; close all; clc

%% Create data vectors

x=[0.1 0.4 0.7 1.2 1.3 1.7 2.2 2.8 3.0 4.0 4.3 4.4 4.9];
y=[0.5 0.9 1.1 1.5 1.5 2.0 2.2 2.8 2.7 3.0 3.5 3.7 3.9];

%% Best fit by error function

line = @ (m,b) m*x+b;
E1 = @ (par) norm(y - line(par(1),par(2)),1); % par = parameter = (m,b) 
E2 = @ (par) norm(y - line(par(1),par(2)),2);
E10 = @ (par) norm(y - line(par(1),par(2)),10); % just for fun
coeff_L1 = fminsearch(E1, [1 1]); % Initial guess for par = [1 1]
coeff_L2 = fminsearch(E2, [1 1]);
coeff_L10 = fminsearch(E10, [1 1]);

%% Plot results

figure(1)
plot(x,y,'ok','Linewidth',2)
hold on
plot(x,line(coeff_L1(1),coeff_L1(2)),'b','Linewidth',2)
plot(x,line(coeff_L2(1),coeff_L2(2)),'r','Linewidth',2)
plot(x,line(coeff_L10(1),coeff_L10(2)),'g','Linewidth',2)
legend('data','L^1 fit','L^2 fit', 'L^{10} fit', 'Location', 'southeast')
set(gca,'Fontsize',16)

%% Add in outliers

x = [x 0.5 0.8];
y = [y 3.9 0.3];

line = @ (m,b) m*x+b;
E1 = @ (par) norm(y - line(par(1),par(2)),1); % par = parameter = (m,b) 
E2 = @ (par) norm(y - line(par(1),par(2)),2);
E10 = @ (par) norm(y - line(par(1),par(2)),10); % just for fun
coeff_L1 = fminsearch(E1, [1 1]); % Initial guess for par = [1 1]
coeff_L2 = fminsearch(E2, [1 1]);
coeff_L10 = fminsearch(E10, [1 1]);

figure(2)
plot(x,y,'ok','Linewidth',2)
hold on
plot(x,line(coeff_L1(1),coeff_L1(2)),'b','Linewidth',2)
plot(x,line(coeff_L2(1),coeff_L2(2)),'r','Linewidth',2)
plot(x,line(coeff_L10(1),coeff_L10(2)),'g','Linewidth',2)
legend('data','L^1 fit','L^2 fit', 'L^{10} fit', 'Location', 'southeast')
set(gca,'Fontsize',16)

%% Underdetermined linear system

clear all; close all; clc

A = [1/2 1/4 3 7 -2 10];
b = [3];
x = A\b;

%% Solving with the pseudoinverse

x2 = pinv(A)*b;

%% Comparing solution finding methods

% Initializing a 100x500 linear system
m=100; n=500;
A=randn(m,n);
b=randn(m,1);

% Built-in MATLAB solvers
x1=A\b;
x2=pinv(A)*b;

% cvx solutions
cvx_begin quiet
    variable x3(n);
    minimize( norm(x3,1) );
    subject to
    A*x3 == b;
cvx_end
cvx_begin quiet
    variable x4(n);
    minimize( norm(x4,2) );
    subject to
    A*x4 == b;
cvx_end

% Plot histograms
figure(3)
subplot(4,1,1)
histogram(x1,linspace(-.2,.2,31)) % specify edges of bins
title('Backslash')
xlim([-.2 .2])
subplot(4,1,2)
histogram(x2,linspace(-.2,.2,31))
title('Pseudoinverse')
xlim([-.2 .2])
subplot(4,1,3)
histogram(x1,linspace(-.2,.2,31))
title('L1 minimization')
xlim([-.2 .2])
subplot(4,1,4)
histogram(x4,linspace(-.2,.2,31))
title('L2 minimization')

%% Number of nonzero elements in our solutions

nnz(x1)
nnz(x3)

%% Solving an overdetermined system

clear all; close all; clc

m=500; n=150;
A=randn(m,n);
b=randn(m,1);

% Built-in MATLAB solvers
x1=A\b;
x2=pinv(A)*b;

% cvx solutions
cvx_begin quiet
    variable x3(n);
    minimize( norm(A*x3-b,1) );
cvx_end
cvx_begin quiet
    variable x4(n);
    minimize( norm(A*x4-b,2) );
cvx_end

% Plot historgrams
figure(4)
subplot(4,1,1)
histogram(A*x1-b,80) % specify number of bins
title('Backslash')
xlim([-4 4])
subplot(4,1,2)
histogram(A*x2-b,80)
title('Pseudoinverse')
xlim([-4 4])
subplot(4,1,3)
histogram(A*x3-b,80)
title('L1 minimization')
xlim([-4 4])
subplot(4,1,4)
histogram(A*x4-b,80)
title('L2 minimization')
xlim([-4 4])

