%% Lecture 17 Code

% Clean workspace
clear all; close all; clc

%% Load in sound clips

% Sound files to load
files = {'chirp.mat'
        'gong.mat'
        'handel.mat'
        'laughter.mat'
        'splat.mat'
        'train.mat'};

% Matrix of principal components
S = zeros(10000,6);
for i = 1:6
    test     = load(files{i});
    y        = test.y(1:10000,1);
    S(:,i)   = y;
end

% Play sound
soundsc(S(:,1))

%% Mix the signals together

rng default % For reproducibility
mixdata = S*randn(6) + randn(1,6);

% listen to mixed signal
soundsc(mixdata(:,1))

%% Plot the mixed signals

figure(1)
for i = 1:6
    subplot(2,6,i)
    plot(S(:,i))
    title(['Sound ',num2str(i)])
    set(gca,'Fontsize',12)
    subplot(2,6,i+6)
    plot(mixdata(:,i))
    title(['Mix ',num2str(i)])
    set(gca,'Fontsize',12)
end

%% Whiten the data

% Compute SVD of covariance matrix
[U,Sig] = svd(cov(mixdata));
Sig     = diag(Sig);
Sig     = Sig(:)';

% Compute whitened data
mu = mean(mixdata,1);
mixdata = bsxfun(@minus,mixdata,mu);
mixdata = bsxfun(@times,mixdata*U,1./sqrt(Sig));

%% Perform RICA

q = 6; % number of independent components
Mdl = rica(mixdata,q); 
unmixed = transform(Mdl,mixdata);

%% Plot un-mixed signals

figure(2)
for i = 1:6
    subplot(2,6,i)
    plot(S(:,i))
    title(['Sound ',num2str(i)])
    set(gca,'Fontsize',12)
    subplot(2,6,i+6)
    plot(unmixed(:,i))
    title(['Unmix ',num2str(i)])
    set(gca,'Fontsize',12)
end


%% Reorder and rescale

unmixed = unmixed(:,[2,5,3,6,3,1]);
for i = 1:6
    unmixed(:,i) = unmixed(:,i)/norm(unmixed(:,i))*norm(S(:,i));
end

figure(3)
for i = 1:6
    subplot(2,6,i)
    plot(S(:,i))
    ylim([-1,1])
    title(['Sound ',num2str(i)])
    set(gca,'Fontsize',12)
    subplot(2,6,i+6)
    plot(unmixed(:,i))
    ylim([-1,1])
    title(['Unmix ',num2str(i)])
    set(gca,'Fontsize',12)
end






