function [deno, aggw, bias] = nlbayes_step(nisy, bsic, orig, sigma, prms)


% binary that computes patch distances
pgbin = '/home/pariasm/Work/denoising/projects/video_nlbayes3d/build/bin/patch_group';

% save sequence
ori_pat = '/tmp/tmp_ori_%03d.png';
bsc_pat = '/tmp/tmp_bsc_%03d.png';

% in the absence of a bias
if isempty(bsic),
	bsic = nisy;
	step2 = false;
else
	step2 = true;
end

% seqs are saved to input them to the c++ distance computation
save_sequence(bsic, bsc_pat, 1);
if ~isempty(orig),
	save_sequence(orig, ori_pat, 1);
end

% command to compute nearest patches
command = [pgbin ...
           ' -i ' bsc_pat ...
           ' -f ' num2str(1) ...
           ' -l ' num2str(size(bsic,4)) ...
           ' -sigma ' num2str(sigma) ...
           ' -px ' num2str(prms.px) ...
           ' -pt ' num2str(prms.pt) ...
           ' -wx ' num2str(prms.wx) ...
           ' -wt ' num2str(prms.wt) ...
           ' -np ' num2str(prms.np)   ];

%command = [pgbin ...
%           ' -i ' ori_pat ...
%           ' -f ' num2str(1) ...
%           ' -l ' num2str(size(orig,4)) ...
%           ' -sigma ' sigma ...
%           ' -px ' num2str(prms.px) ...
%           ' -pt ' num2str(prms.pt) ...
%           ' -wx ' num2str(prms.wx) ...
%           ' -wt ' num2str(prms.wt) ...
%           ' -np ' num2str(prms.np)   ];


chnls = size(nisy,3);

aggw = zeros(size(nisy));
aggp = zeros(size(nisy));
aggu = zeros(size(nisy));
deno = zeros(size(nisy));

bias = [];
if ~isempty(orig),
	aggb = zeros(size(nisy));
	bias = zeros(size(nisy));
end

stepx = 1;%floor(prms.px/2);
stepy = 1;%floor(prms.px/2);
stept = 1;
ii = 0;

U = [];
S = [];


% build DCT basis
cx = cos((1/2 + [0:prms.px-1]') * [0:prms.px-1]*pi/prms.px);
ct = cos((1/2 + [0:prms.pt-1]') * [0:prms.pt-1]*pi/prms.pt);

cx(:,1) = cx(:,1  ) / sqrt(prms.px); 
ct(:,1) = ct(:,1  ) / sqrt(prms.pt); 
cx(:,2:prms.px) = cx(:,2:prms.px) * sqrt(2/prms.px); 
ct(:,2:prms.pt) = ct(:,2:prms.pt) * sqrt(2/prms.pt); 

U = zeros(prms.px*prms.px*prms.pt, prms.px*prms.px*prms.pt);


for k = 1:prms.pt,
for i = 1:prms.px,
for j = 1:prms.px,

	u = cx(:,i)*cx(:,j)';
	u = u(:)*ct(:,k)';
	U(:,(k-1)*prms.px*prms.px + (i-1)*prms.px + j) = u(:);

end
end
end

for pat = 1:stept:size(nisy,4) - prms.pt+1,
for pay = 1:stepy:size(nisy,1) - prms.px+1,
for pax = 1:stepx:size(nisy,2) - prms.px+1, ii = ii + 1;

	if aggp(pay,pax,1,pat),
		continue;
	end

	% -------------------------------------------------- compute patch group
	cmd = [command ...
	      ' -PAx ' num2str(pax - 1) ... 
	      ' -PAy ' num2str(pay - 1) ... 
	      ' -PAt ' num2str(pat - 1) ... 
	      ' -out /tmp/patch_group.out'];

	unix('rm /tmp/patch_group.out');
	unix(cmd);

	% separate command outputs into matlab readable files
	unix('cat /tmp/patch_group.out | grep COORD | sed -s ''s/\[COORD\]//'' > /tmp/patches.coords.out');
%	unix('cat /tmp/patch_group.out | grep DISTA | sed -s ''s/\[DISTA\]//'' > /tmp/patches.distas.out');
%	unix('cat /tmp/patch_group.out | grep YUV0  | sed -s ''s/\[YUV0\]//''  > /tmp/patches.y.out');
%	unix('cat /tmp/patch_group.out | grep YUV1  | sed -s ''s/\[YUV1\]//''  > /tmp/patches.u.out');
%	unix('cat /tmp/patch_group.out | grep YUV2  | sed -s ''s/\[YUV2\]//''  > /tmp/patches.v.out');
%	unix('cat /tmp/patch_group.out | grep BSIC  | sed -s ''s/\[BSIC\]//''  > /tmp/patches.bsic.out');
%	unix('cat /tmp/patch_group.out | grep NISY  | sed -s ''s/\[NISY\]//''  > /tmp/patches.nisy.out');

	patches.coords = load('/tmp/patches.coords.out') + 1;
%	patches.distas = load('/tmp/patches.distas.out');

	group_sz = size(patches.coords,1);

	% ------------------------------------------------------ extract patches
	patches.nisy = zeros(prms.px, prms.px, chnls, prms.pt, group_sz);
	for i = 1:group_sz,
		patches.nisy(:,:,:,:,i) = nisy(patches.coords(i,2):patches.coords(i,2)+prms.px-1,...
		                               patches.coords(i,1):patches.coords(i,1)+prms.px-1,:,...
		                               patches.coords(i,3):patches.coords(i,3)+prms.pt-1);

		if step2,
			patches.bsic(:,:,:,:,i) = bsic(patches.coords(i,2):patches.coords(i,2)+prms.px-1,...
													 patches.coords(i,1):patches.coords(i,1)+prms.px-1,:,...
													 patches.coords(i,3):patches.coords(i,3)+prms.pt-1);
		end

	end

	% bias: extract patches fromo original video
	if ~isempty(orig),
		patches.orig = zeros(prms.px, prms.px, chnls, prms.pt, group_sz);
		for i = 1:group_sz,
			patches.orig(:,:,:,:,i) = orig(patches.coords(i,2):patches.coords(i,2)+prms.px-1,...
													 patches.coords(i,1):patches.coords(i,1)+prms.px-1,:,...
													 patches.coords(i,3):patches.coords(i,3)+prms.pt-1);
		end
	else
		patches.orig = [];
	end

	% ----------------------------------------------- compute bayes estimate
	gsz  = size(patches.nisy);
	pdim = prod(gsz(1:4));

	nn = reshape(patches.nisy,[pdim gsz(5)]);
	oo = [];
	if ~isempty(orig),
		oo = reshape(patches.orig,[pdim gsz(5)]);
	end

	if step2,
		aa = reshape(patches.bsic,[pdim gsz(5)]);
%		U = [];
		S = [];
		[dd,m,U,S,bb] = compute_bayes_estimate(nn,aa,sigma,prms.r,oo,'pos', U,S);
	else
%		U = [];
		S = [];
		[dd,m,U,S,bb] = compute_bayes_estimate(nn,[],sigma,prms.r,oo,prms.filter_type, U,S);
	end

	dd = reshape(dd,gsz);
	if ~isempty(orig),
		bb = reshape(bb,gsz);
	end


	% ------------------------------------------- aggregate patches on image
	for i = 1:group_sz,

		cc = patches.coords;

%		if aggp(cc(i,2), cc(i,1), 1, cc(i,3)) == 0, 

		aggu(cc(i,2):cc(i,2)+prms.px-1,...
		     cc(i,1):cc(i,1)+prms.px-1,:,...
		     cc(i,3):cc(i,3)+prms.pt-1) = dd(:,:,:,:,i) + ...
		                             1*aggu(cc(i,2):cc(i,2)+prms.px-1,...
		                                  cc(i,1):cc(i,1)+prms.px-1,:,...
		                                  cc(i,3):cc(i,3)+prms.pt-1); 

		aggw(cc(i,2):cc(i,2)+prms.px-1,...
		     cc(i,1):cc(i,1)+prms.px-1,:,...
		     cc(i,3):cc(i,3)+prms.pt-1) = 1 + ...
		                             1*aggw(cc(i,2):cc(i,2)+prms.px-1,...
		                                  cc(i,1):cc(i,1)+prms.px-1,:,...
		                                  cc(i,3):cc(i,3)+prms.pt-1); 

		if ~isempty(orig),
			aggb(cc(i,2):cc(i,2)+prms.px-1,...
			     cc(i,1):cc(i,1)+prms.px-1,:,...
			     cc(i,3):cc(i,3)+prms.pt-1) = bb(:,:,:,:,i) + ...
			                             1*aggb(cc(i,2):cc(i,2)+prms.px-1,...
			                                  cc(i,1):cc(i,1)+prms.px-1,:,...
			                                  cc(i,3):cc(i,3)+prms.pt-1); 
		end

		aggp(patches.coords(i,2),...
		     patches.coords(i,1),:,...
		     patches.coords(i,3)) = 1 + ...
		                          1*aggp(patches.coords(i,2),...
		                               patches.coords(i,1),:,...
		                               patches.coords(i,3)); 

%		end

%		nonzero = find(aggw ~= 0);
%		deno(nonzero) = min(255, max(0, aggu(nonzero) ./ aggw(nonzero)));
		

		%figure(1)
%		imagesc(deno(:,:,:,pat+floor(prms.pt/2))/255);
%		axis equal, axis off,
%		drawnow
%		pause(.01)

	end

	if (mod(ii,5) == 1)
		nonzero = find(aggw ~= 0);
		deno_prev = deno;
		deno(nonzero) = min(255, max(0, aggu(nonzero) ./ aggw(nonzero)));

		if ~isempty(orig),
			bias(nonzero) = min(255, max(0, aggb(nonzero) ./ aggw(nonzero)));
		end

%		disp(sqrt(norm(deno(:) - deno_prev(:))^2/prms.np/prms.px/prms.pt/chnls))

	%	figure(1)

		pat_show = pat+floor(prms.pt/2);
		if isempty(orig),
			imagesc([deno(:,:,:,pat_show)/255 aggw(:,:,:,pat_show)/pdim/prms.np*10],[0 1]);
		else
			imagesc([bias(:,:,:,pat_show)/255 deno(:,:,:,pat_show)/255],[0 1]);
		end
		axis equal, axis off,
		colormap gray
		drawnow
		pause(.01)
	end

end
end
end



