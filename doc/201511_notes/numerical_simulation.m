% This script generates follows the notes. We generate M random observations
% of N vectors x. These N vectors are generated from a normal diststribution
% with diagonal covariance matrix. 
%
% For each of the M batches, we compute MAP estimates of the x vectors using
% an empirical Wiener filter. These M estimates are collected, so that we
% can compute their empirical mean and variance, and check if they correspond
% to the mean and variance we derived in the notes.

% ---------------------------------------------------------- create model

d = 100;						% # of observed variables (dimension of ambient space)
N = 200;						% # of data points
M = 100;						% # of aggregated estimates

% Square root of covaciance matrix
W  = 2*diag([d-1:-1:0]);

% Noise
sigma = 20;

% ---------------------------------------------- generate data from model
Z = randn(d,N); 
X0 = W*Z;


% ------------------------------------------------------ aggregation loop

X_map1 = 0;
X_map2 = 0;
X_map1s = zeros(d,N,M);
X_map2s = zeros(d,N,M);

X_mean = zeros(d,N);

for aa = 1:M,

	noise = sigma*randn(d,N);
	X = X0 + noise;

	% ------------------------------------------------------ learn parameters
	
	% maximum likelihood estimates for noise (known) and mean
	S = 1/N*sum(X.^2,2);
	
	% ---------------------------------------------------------- estimate MAP
	sigma_agg = sigma;

	% estimate MAPS are stored
	X_map1s(:,:,aa) = diag(    1 - sigma_agg^2./S   ) * X;
	X_map2s(:,:,aa) = diag(max(1 - sigma_agg^2./S,0)) * X;
	
	% aggregate MAPs estimate
	X_map1 = (aa - 1)*X_map1 + X_map1s(:,:,aa); X_map1 = X_map1 / aa;
	X_map2 = (aa - 1)*X_map2 + X_map2s(:,:,aa); X_map2 = X_map2 / aa;

	% compute sample mean
	X_mean = (aa - 1)*X_mean + X; X_mean = X_mean / aa;
	
%	if mod(aa,5) == 0, 
%		figure(2), clf, hold on
%		plot([1:d], X0(:,1), '.-',...
%			  [1:d], X_map1(:,1), '.-',...
%			  [1:d], X_map2(:,1), '.-')
%		legend('gt','full rank w/o th', 'full rank w th')
%		
%		pause(0.1)
%	end
	
end

% empirical variance of the MAP estimates
V_map1 = mean((X_map1s - repmat(X_map1,[1 1 M])).^2,3);
V_map2 = mean((X_map2s - repmat(X_map2,[1 1 M])).^2,3);

% empirical bias
b_map1 = mean((X_map1s - repmat(X0    ,[1 1 M])).^2,3) - V_map1;
b_map2 = mean((X_map2s - repmat(X0    ,[1 1 M])).^2,3) - V_map2;

% empirical signal-to-noise ratio
snr = mean(X0.^2,2)/sigma^2;

% approx expected value of Wiener coefficients and MAP estimates
E_w = 1./(1 + 1./snr);
E_x = diag(E_w)*X0;

% approx variance of Wiener coefficients and MAP estimates
V_w = 2/N*(1 + 2*snr)./(1 + snr).^4;
V_x = diag(V_w)*(X0.^2 + sigma^2) + sigma^2*E_w.^2*ones(1,N);

% true snr and corresponding theoretical variance and bias
snr_t = diag(W).^2/sigma^2;
V_x_t = 2/N*sigma^2*(1 + 2*snr_t)./(1 + snr_t).^3 + sigma^2*snr_t.^2./(1 + snr_t).^2;
b_t   = sigma^2 * snr_t./(snr_t + 1).^2;


% bias of MAP estimates
b = diag(1./(1 + snr))*X0;

% plot: compare the estimation of a vector of the set
figure(1)
plot(1:d, X0(:,1),'.-', 1:d,E_x(:,1),'.-', 1:d, X_map1(:,1),'.-', 1:d, X_map2(:,1), 1:d, X_mean(:,1),'.-')
h = legend('$x_i$', '$E[\hat x_i]$', '$\hat{x}_i$', '$\hat x_i^+$', 'Sample average');
set(h,'Interpreter','latex');
title('Aggregated estimate and expected value of sigle estimate','Interpreter','latex')
xlabel('$k$','Interpreter','latex')

% plots: theoretical and empirical variances, bias and sample mean variance
figure(2)
%plot(1:d,V_x(:,1),'.-', 1:d, V_map1(:,1),'.-', 1:d, V_map2(:,1),'.-', ...
%     [1 d], [1 1]*sigma^2, 'k--' , 1:d, b(:,1).^2, '.-', 1:d, V_x_t, '--', 1:d, b_t)
%h = legend('$MV[\hat x_i]$', '$M\widehat{V}[\hat x_i]$', '$M\widehat{V}[\hat x_i^+]$',...
%           '$\sigma^2$', '$1/N \sum_{i=1}^N B[x_i]^2$','Location','NorthWest');
plot(1:d,mean(V_map1,2),'.-', 1:d,mean(V_map2,2),'.-', ...
     1:d,mean(b_map1,2),'.-', 1:d,mean(b_map2,2),'.-', ...
     [1 d], [1 1]*sigma^2, 'k--' ,...
	  1:d, V_x_t, '--', 1:d, b_t)
ylim([0 sigma^2 + 50])
grid on
h = legend('$\frac MN \sum_{i=1}^N \widehat V\{\hat x_i\}$', ...
           '$\frac MN \sum_{i=1}^N \widehat V\{\hat x^+_i\}$', ...
           '$\frac 1N \sum_{i=1}^N \widehat B\{\hat x_i\}$', ...
           '$\frac 1N \sum_{i=1}^N \widehat B\{\hat x^+_i\}$', ...
           '$\sigma^2$',...
			  '$\frac MN \sum_{i=1}^N V\{\hat x_i\}$',...
			  '$\frac 1N \sum_{i=1}^N B\{\hat x_i\}$',...
			  'Location','West');
set(h,'Interpreter','latex');
%title('Variance of the single estimates','Interpreter','latex')
xlabel('$k$','Interpreter','latex')


disp('RMSE of MAP with coefficient thresholding:')
disp(sqrt(norm(X0 - X_map2, 'fro').^2/d/N))

disp('RMSE of MAP w/o thresholding:')
disp(sqrt(norm(X0 - X_map1, 'fro').^2/d/N))

disp('RMSE of MAP w/o thresholding using theoretical bias-variance:')
disp(sqrt(mean(mean(b.^2,2) + mean(V_x,2)/M)))

disp('RMSE of sample average:')
disp(sqrt(norm(X0 - X_mean, 'fro').^2/d/N))

