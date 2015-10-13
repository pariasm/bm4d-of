function [deno, mean_nisy, U, S, bias] = compute_bayes_estimate(nisy, bsic, sigma, rank, orig, filter_type, U, S)

if nargin < 6,
	filter_type = 'pos';
end

if nargin < 7,
	U = [];
	S = [];
end

d = size(nisy,1);
n = size(nisy,2);

sigma_bsic = 0;
if isempty(bsic),
	bsic = nisy;
	sigma_bsic = sigma;
end


% center
%mean_nisy = 128 + 0*mean(nisy, 2);
%mean_bsic = 128 + 0*mean(bsic, 2);
mean_nisy = mean(nisy, 2);
mean_bsic = mean(bsic, 2);
nisy = nisy - mean_nisy*ones(1,n);
bsic = bsic - mean_bsic*ones(1,n);

if ~isempty(orig),
	orig = orig - mean_nisy*ones(1,n);
end


if isempty(S),

	if isempty(U),
		% covariance matrix
		C = 1/n*bsic*bsic';
		
		% eigendecomposition
		if rank > 0 && rank < d,
			[U,S] = eigs(C,rank);
		else
			[U,S] = eig(C);
		end
		
		U = real(U);
		S = real(S);
	else
		% project data vectors onto basis
		S = diag(mean((U'*bsic).^2,2));
	end
	
	
	% sort
	tmp = sortrows([ diag(S) U'],-1);
	S = diag(tmp(:,1));
	U = tmp(:,2:end)';
	
	% wiener filter coefficients
	%S = max(0, diag(S) - sigma_bsic*sigma_bsic);
	if     strcmp(filter_type, 'pos'), 
		S = min(max( 0  , diag(S) - sigma_bsic*sigma_bsic),Inf);
		S = diag(1./(1 + sigma*sigma ./ S));
%		S = diag(double(S > 0));
%		S = 0*diag(double(S > 0));
	elseif strcmp(filter_type, 'neg'),
		S = min(max(-Inf, diag(S) - sigma_bsic*sigma_bsic),Inf);
		S = diag(1./(1 + sigma*sigma ./ S));
%		S = diag(2*double(S > 0) - 1);
%		S = 0*diag(double(S > 0));
	elseif strcmp(filter_type, 'neg-inv'),
		S = min(max(-Inf, diag(S) - sigma_bsic*sigma_bsic),Inf);
		S = abs(diag(1./(1 + sigma*sigma ./ S)));
	end

elseif isempty(S)




end

%i = 196;

% denoise group
deno = (U*S)*(U'*nisy) + mean_nisy*ones(1,n);
%deno = (U*S)*(U'*nisy) + 128;

bias = [];
if ~isempty(orig),
	bias = (U*S)*(U'*orig) + mean_nisy*ones(1,n);
%	bias = (U*S)*(U'*orig) + 128;
end
%deno = 5*(U*S)*(U'*nisy) + 128;
%deno = (U(:,i)*S(i,i))*(U(:,i)'*nisy) + mean_nisy*ones(1,n);
%deno = (U(:,i)*i*S(i,i))*(U(:,i)'*nisy) + 128;



