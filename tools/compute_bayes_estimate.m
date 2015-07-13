function [deno, mean_bsic, U, S] = compute_bayes_estimate(nisy, bsic, sigma, rank)

d = size(nisy,1);
n = size(nisy,2);

sigma_bsic = 0;
if isempty(bsic),
	bsic = nisy;
	sigma_bsic = sigma;
end


% center
mean_nisy = mean(nisy, 2);
mean_bsic = mean(bsic, 2);
nisy = nisy - mean_nisy*ones(1,n);
bsic = bsic - mean_bsic*ones(1,n);


% covariance matrix
C = 1/n*bsic*bsic';


% eigendecomposition
[U,S] = eigs(C,rank);


% sort well
tmp = sortrows([ diag(S) U'],-1);
S = diag(tmp(:,1));
U = tmp(:,2:end)';


% wiener filter coefficients
S = max(0, diag(S) - sigma_bsic*sigma_bsic);
S = diag(1./(1 + sigma*sigma ./ S));


% denoise group
deno = (U*S)*(U'*nisy) + mean_nisy*ones(1,n);


S = diag(S);


