% loads a video and allows the user to select a pixel
% for the selected pixel, it displays the locations 
% of the most similar patches

%clear all


% parameters
prms.wx = 201;
prms.wt = 0;
prms.px = 4;
prms.pt = 1;
prms.np = 200;
sigma = 20;

% path to data
%base_path = '../../../data/';
base_path = '../build/';

% image data
%orig_path = [base_path 'images/tune/'];
%seq.name = 'Girl.png';
orig_path = base_path;
seq.name = 'bsic41.png';
seq.first = 0;
seq.last  = 0;
%seq.pattern = [orig_path 'girl.png'];
seq.pattern = [orig_path 'bsic41.png'];

% % sequence data
% orig_path = [base_path 'data/derf/'];
% seq.name = 'bus';
% seq.first = 1;
% seq.last  = 150;
% seq.pattern = [orig_path seq.name '/%03d.png'];

% load
orig = load_sequence(seq.pattern, seq.first, seq.last);
nisy = orig + 0*sigma*randn(size(orig));
nisy0 = nisy;
orig0 = orig;

[HH WW CC FF] = size(nisy);

% show sequence
goon = 1;
frame = 1;

main_fig_h = figure(1);
imagesc(uint8(nisy(:,:,:,frame)));
main_axes_h = gca;
axis equal, axis off

key = '';
set(main_fig_h,'KeyPressFcn', @(~, evt) assignin('base', 'key', evt.Key)); 


% build DCT basis
cx = cos((1/2 + [0:prms.px-1]') * [0:prms.px-1]*pi/prms.px);
ct = cos((1/2 + [0:prms.pt-1]') * [0:prms.pt-1]*pi/prms.pt);

cx(:,1) = cx(:,1  ) / sqrt(prms.px); 
ct(:,1) = ct(:,1  ) / sqrt(prms.pt); 
cx(:,2:prms.px) = cx(:,2:prms.px) * sqrt(2/prms.px); 
ct(:,2:prms.pt) = ct(:,2:prms.pt) * sqrt(2/prms.pt); 

U = zeros(prms.px*prms.px*prms.pt, prms.px*prms.px*prms.pt);

% window
wx = chebwin(prms.px);
wwx = repmat(reshape(wx*wx',[prms.px*prms.px 1]),[3 1]);


for k = 1:prms.pt,
for i = 1:prms.px,
for j = 1:prms.px,

	u = cx(:,i)*cx(:,j)';
	u = u(:)*ct(:,k)';
	U(:,(k-1)*prms.px*prms.px + (i-1)*prms.px + j) = u(:);

end
end
end

U = kron(eye(CC), U);

ind = [];
while (goon)

	disp(key)

	switch (key),
	case 'q'
		goon = 0;

	case 'f',

		% advance to next frame
		frame = min(frame + 1, size(nisy,4));

		figure(main_fig_h)
		current_axes = [get(gca,'XLim') get(gca,'YLim')];
		imagesc(uint8(nisy(:,:,:,frame)),'parent',main_axes_h);
		axis equal, axis off, axis(current_axes)
		hold on
		%plot(pax, pay, '*r')
		hold off

	case 'b',

		% backward to previous frame
		frame = max(frame - 1, 1);

		figure(main_fig_h)
		current_axes = [get(gca,'XLim') get(gca,'YLim')];
		imagesc(uint8(nisy(:,:,:,frame)),'parent',main_axes_h);
		axis equal, axis off, axis(current_axes)
		hold on
		%plot(pax, pay, '*r')
		hold off

	case {'s','c'},

		if (key == 's'),
			% select a pixel
			figure(main_fig_h)
			[pax,pay] = ginput(1);
			pax = round(pax);
			pay = round(pay);
			pat = frame;
		end

		% compute patch group
		srch_region = nisy(max(1,pay - prms.wx):min(HH,pay + prms.wx + prms.px - 1),...
								 max(1,pax - prms.wx):min(WW,pax + prms.wx + prms.px - 1),:,...
								 max(1,pat - prms.wt):min(FF,pat + prms.wt + prms.pt - 1));

		srch_patches = vid2col_ch(srch_region, [prms.px prms.px prms.pt]);
		refe_patch = nisy(pay:pay + prms.px - 1,...
		                  pax:pax + prms.px - 1,:,...
		                  pat:pat + prms.pt - 1);

		[patches.distas, idx] = sort(L2_distance(refe_patch(:), srch_patches));

		% coordinates of the np nearest neighbors to the ref patch
		idx = idx(1:prms.np)';

		srch_h  =  size(srch_region,1) - prms.px + 1;
		srch_wh = (size(srch_region,2) - prms.px + 1)*srch_h;
		patches.coords = [max(1,pax - prms.wx) + floor(mod(idx-1, srch_wh)/srch_h),...
								max(1,pay - prms.wx) +   mod(    idx-1,          srch_h),...
								max(1,pat - prms.wt) + floor(   (idx-1)/srch_wh)       ];


		% extract similar patches
		patches.nisy = srch_patches(:,idx);

		srch_region = orig(max(1,pay - prms.wx):min(HH,pay + prms.wx + prms.px - 1),...
								 max(1,pax - prms.wx):min(WW,pax + prms.wx + prms.px - 1),:,...
								 max(1,pat - prms.wt):min(FF,pat + prms.wt + prms.pt - 1));
		patches.orig = vid2col_ch(srch_region, [prms.px prms.px prms.pt]);
		patches.orig = patches.orig(:,idx);


		% on a separate figure, display set of similar patches
		np_viz = 20;

		npatches = size(patches.nisy,2);
		patchdim = size(patches.nisy,1);
		rank = 40;
		[dd,m] = compute_bayes_estimate(patches.nisy,[],sigma,rank,[],'neg',U,[],wwx);

		dd = min(255, max(0, dd));
		oo = patches.orig;
		pp = patches.nisy;

		pp = build_patch_image(reshape(pp(:,1:np_viz),[prms.px prms.px CC prms.pt np_viz]));
		dd = build_patch_image(reshape(dd(:,1:np_viz),[prms.px prms.px CC prms.pt np_viz]));
		oo = build_patch_image(reshape(oo(:,1:np_viz),[prms.px prms.px CC prms.pt np_viz]));

		figure(2)
		imagesc(uint8([pp ; dd ; oo])), axis equal, axis off

		% display also eigenvectors and eigenvalues
		mU = min(255,max(0,[m, 255*(4*U(:,1:np_viz-1) + 0.5)]));
		mU = build_patch_image(reshape(mU,[prms.px prms.px 3 prms.pt np_viz]));
		figure(3)
		imagesc(mU/255), axis equal, axis off




		% show search region with selected patches
		indr = sub2ind(size(nisy),patches.coords(:,2), ...
		                          patches.coords(:,1), ...
		                          ones(size(patches.coords,1),1), ...
		                          patches.coords(:,3));

		indg = sub2ind(size(nisy),patches.coords(:,2), ...
		                          patches.coords(:,1), ...
		                          2*ones(size(patches.coords,1),1), ...
		                          patches.coords(:,3));

		indb = sub2ind(size(nisy),patches.coords(:,2), ...
		                          patches.coords(:,1), ...
		                          3*ones(size(patches.coords,1),1), ...
		                          patches.coords(:,3));

		a = 0.0;
		tmp = 0.5*repmat(mean(nisy,3),[1 1 3]) + 0.5*nisy;
		tmp(indr(1 :  5)) = a*tmp(indr(1 : 5 )) + (1-a)*255;
		tmp(indr(6 :end)) = a*tmp(indr(6 :end)) + (1-a)*  0;
%		tmp(indr(6 : 45)) = a*tmp(indr(6 :45 )) + (1-a)*  0;
%		tmp(indr(46:end)) = a*tmp(indr(46:end)) + (1-a)*  0;
		tmp(indg(1 : 5 )) = a*tmp(indg(1 : 5 )) + (1-a)*  0;
		tmp(indg(6 :end)) = a*tmp(indg(6 :end)) + (1-a)*255;
%		tmp(indg(6 :45 )) = a*tmp(indg(6 :45 )) + (1-a)*255;
%		tmp(indg(46:end)) = a*tmp(indg(46:end)) + (1-a)*  0;
		tmp(indb(1 : 5 )) = a*tmp(indb(1 : 5 )) + (1-a)*  0;
		tmp(indb(6 :end)) = a*tmp(indb(6 :end)) + (1-a)*  0;
%		tmp(indb(6 :45 )) = a*tmp(indb(6 :45 )) + (1-a)*  0;
%		tmp(indb(46:end)) = a*tmp(indb(46:end)) + (1-a)*255;
		search.rx = floor(prms.wx/2);
		search.rt =       prms.wt;
		search.box_x = [max(1, pax - search.rx - prms.px), min(size(nisy,2), pax + search.rx + prms.px)];
		search.box_y = [max(1, pay - search.rx - prms.px), min(size(nisy,1), pay + search.rx + prms.px)];
		search.box_t = [max(1, pat - search.rt          ), min(size(nisy,4), pat + search.rt          )];
		search.image = tmp(search.box_y(1): search.box_y(2),...
		                   search.box_x(1): search.box_x(2),:,...
		                   search.box_t(1): search.box_t(2));

		aa = search.image;
		sza = size(aa);
		if (length(sza) == 3) sza(4) = 1; end
		aa = reshape( cat(2, permute(aa, [1 2 4 3]), 255*ones(sza(1),3, sza(4), sza(3))),...
		              [sza(1), sza(4)*(sza(2) + 3), sza(3)]);
		aa = aa(:,1:end-3,:);

		figure(4)
		imagesc(uint8(aa))
		axis equal, axis off

%		% save results
%		fname = sprintf('patch_group_%s_%03d_%03d_%03d_s%02d_wx%d_wt%d_sx%d_st%d_r%03d_n%03d_pcas.png',seq.name, pax,pay,pat, sigma, prms.wx,prms.wt,prms.px,prms.pt, rank, prms.np);
%		imwrite(uint8(mU),fname)                                                                                            
%		fname = sprintf('patch_group_%s_%03d_%03d_%03d_s%02d_wx%d_wt%d_sx%d_st%d_r%03d_n%03d_coor.png',seq.name, pax,pay,pat, sigma, prms.wx,prms.wt,prms.px,prms.pt, rank, prms.np);
%		imwrite(uint8(aa/255*255),fname)                                                                                            
%		fname = sprintf('patch_group_%s_%03d_%03d_%03d_s%02d_wx%d_wt%d_sx%d_st%d_r%03d_n%03d_orig.png',seq.name, pax,pay,pat, sigma, prms.wx,prms.wt,prms.px,prms.pt, rank, prms.np);
%		imwrite(uint8(oo),fname)                                                                                            
%		fname = sprintf('patch_group_%s_%03d_%03d_%03d_s%02d_wx%d_wt%d_sx%d_st%d_r%03d_n%03d_nisy.png',seq.name, pax,pay,pat, sigma, prms.wx,prms.wt,prms.px,prms.pt, rank, prms.np);
%		imwrite(uint8(pp),fname)                                                                                            
%		fname = sprintf('patch_group_%s_%03d_%03d_%03d_s%02d_wx%d_wt%d_sx%d_st%d_r%03d_n%03d_deno.png',seq.name, pax,pay,pat, sigma, prms.wx,prms.wt,prms.px,prms.pt, rank, prms.np);
%		imwrite(uint8(dd),fname)



		figure(main_fig_h)
		current_axes = [get(gca,'XLim') get(gca,'YLim')];
		imagesc(uint8(nisy(:,:,:,frame)),'parent',main_axes_h);
		axis equal, axis off, axis(current_axes)
		hold on
		%plot(pax, pay, '*r')
		hold off


		% restore original values
		nisy(ind) = nisy0(ind);

	otherwise,
	end


	key = '';

	pause(.1)
end



