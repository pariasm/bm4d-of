% loads a video and allows the user to select a pixel
% for the selected pixel, it displays the locations 
% of the most similar patches

%clear all


% parameters
prms.wx = 37;
prms.wt = 2;
prms.px = 9;
prms.pt = 4;
prms.np = 160;
sigma = '0';

% binary that computes patch distances
pgbin = '/home/pariasm/Work/denoising/projects/video_nlbayes3d/build/bin/patch_group';

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37_2r16_2np160/'];
orig_path = [base_path 'data/derf/'];

% sequence data
seqnames = {'bus','coastguard','foreman','tennis'};
first    = {1    , 1          , 1       , 1      };
last     = {150  , 300        , 300     , 150    };

sequences = struct('name',seqnames,'first',first,'last',last);
clear seqnames first last

% select a sequence
seq = sequences(1);

ori_pat = [orig_path seq.name '/%03d.png'];                     % ground truth
if strcat(sigma,'0'),
	nsy_pat = ori_pat; % noisy
else
	nsy_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png']; % noisy
end
bsc_pat = [nlb3_path seq.name '_s' sigma '_pt4/bsic_%03d.png']; % nlb3 pt2

command = [pgbin ...
           ' -i ' ori_pat ...
           ' -b ' ori_pat ...
           ' -f ' num2str(seq.first) ...
           ' -l ' num2str(seq.last) ...
           ' -sigma ' sigma ...
           ' -px ' num2str(prms.px) ...
           ' -pt ' num2str(prms.pt) ...
           ' -wx ' num2str(prms.wx) ...
           ' -wt ' num2str(prms.wt) ...
           ' -np ' num2str(prms.np)   ];

% -PAx 100 -PAy 200 -PAt 3  > patches

% load
nisy = load_sequence(nsy_pat, seq.first, seq.last);
orig = load_sequence(ori_pat, seq.first, seq.last);
nisy0 = nisy;
orig0 = orig;

% show sequence
goon = 1;
frame = 1;
%pax = -1;
%pay = -1;
%pat = -1;

figure(1)
imagesc(nisy(:,:,:,frame)/255);
axis equal, axis off

ind = [];
while (goon)

	key = input('','s');

	switch (key),
	case 'f',
		frame = min(frame + 1, size(nisy,4));
	case 'b',
		frame = max(frame - 1, 1);
	case {'s','c'},
		if (key == 's'),
			figure(1)
			[pax,pay] = ginput(1)
			pax = round(pax);
			pay = round(pay);
			pat = frame;
		end

		% compute patch group
		cmd = [command ...
		      ' -PAx ' num2str(pax - 1) ... 
		      ' -PAy ' num2str(pay - 1) ... 
		      ' -PAt ' num2str(pat - 1) ... 
		      ' -out /tmp/patch_group.out'];

		unix('rm /tmp/patch_group.out');
		unix(cmd);

		% separate command outputs into matlab readable files
		unix('cat /tmp/patch_group.out | grep COORD | sed -s ''s/\[COORD\]//'' > /tmp/patches.coords.out');
		unix('cat /tmp/patch_group.out | grep DISTA | sed -s ''s/\[DISTA\]//'' > /tmp/patches.distas.out');
		unix('cat /tmp/patch_group.out | grep YUV0  | sed -s ''s/\[YUV0\]//''  > /tmp/patches.y.out');
		unix('cat /tmp/patch_group.out | grep YUV1  | sed -s ''s/\[YUV1\]//''  > /tmp/patches.u.out');
		unix('cat /tmp/patch_group.out | grep YUV2  | sed -s ''s/\[YUV2\]//''  > /tmp/patches.v.out');
		unix('cat /tmp/patch_group.out | grep BSIC  | sed -s ''s/\[BSIC\]//''  > /tmp/patches.bsic.out');
		unix('cat /tmp/patch_group.out | grep NISY  | sed -s ''s/\[NISY\]//''  > /tmp/patches.nisy.out');

		patches.coords = load('/tmp/patches.coords.out') + 1;
		patches.distas = load('/tmp/patches.distas.out');

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
		nisy(indr) = 255;

		patches.nisy = zeros(prms.px, prms.px, 3, prms.pt, length(patches.distas));
		for i = 1:length(patches.distas),
			patches.nisy(:,:,:,:,i) = nisy0(patches.coords(i,2):patches.coords(i,2)+prms.px-1,...
			                                patches.coords(i,1):patches.coords(i,1)+prms.px-1,:,...
			                                patches.coords(i,3):patches.coords(i,3)+prms.pt-1);
		end

		patches.orig = zeros(prms.px, prms.px, 3, prms.pt, length(patches.distas));
		for i = 1:length(patches.distas),
			patches.orig(:,:,:,:,i) = orig(patches.coords(i,2):patches.coords(i,2)+prms.px-1,...
			                               patches.coords(i,1):patches.coords(i,1)+prms.px-1,:,...
			                               patches.coords(i,3):patches.coords(i,3)+prms.pt-1);
		end

		% on a separate figure, display set of similar patches
		np_viz = 20;
		pp = patches.nisy(:,:,:,:,1:np_viz);

		gsz  = size(patches.nisy);
		pdim = prod(gsz(1:4));
		[dd,m,U,S] = compute_bayes_estimate(reshape(patches.nisy,[pdim gsz(5)]),[],40,16);
		dd = reshape(dd,gsz);
		dd = min(255, max(0, dd));

		pp = build_patch_image(pp);
		dd = build_patch_image(dd(:,:,:,:,1:np_viz));
		oo = build_patch_image(patches.orig(:,:,:,:,1:np_viz));

		figure(2)
		imagesc([pp ; dd ; oo]/255), axis equal, axis off

		% display also eigenvectors and eigenvalues
		mU = min(255,max(0,[m, 255*(4*U + 0.5)]));
		mU = build_patch_image(reshape(mU,[prms.px prms.px 3 prms.pt 1+size(U,2)]));
		figure(3)
		imagesc(mU/255), axis equal, axis off




		% show search region for all frames
		orig(indr(1: 5))   = 255; orig(indg(1: 5))   =   0; orig(indb(1: 5))   =   0;
		orig(indr(6:20))   =   0; orig(indg(6:20))   = 255; orig(indb(6:20))   =   0;
		orig(indr(21:end)) =   0; orig(indg(21:end)) =   0; orig(indb(21:end)) = 255;
		search.rx = floor(prms.wx/2);
		search.rt =       prms.wt;
		search.box_x = [max(0, pax - search.rx - prms.px), min(size(nisy,2), pax + search.rx + prms.px)];
		search.box_y = [max(0, pay - search.rx - prms.px), min(size(nisy,1), pay + search.rx + prms.px)];
		search.box_t = [max(0, pat - search.rt          ), min(size(nisy,4), pat + search.rt          )];
		search.orig = orig(search.box_y(1): search.box_y(2),...
		                   search.box_x(1): search.box_x(2),:,...
		                   search.box_t(1): search.box_t(2));

		aa = search.orig;
		sza = size(aa);
		aa = reshape( cat(2, permute(aa, [1 2 4 3]), 255*ones(sza(1),3, sza(4), sza(3))),...
		              [sza(1), sza(4)*(sza(2) + 3), sza(3)]);
		aa = aa(:,1:end-3,:);

		figure(4)
		imagesc(aa/255)
		axis equal, axis off

	case 'q'
		goon = 0;
	otherwise,
	end

	figure(1)
	current_axes = [get(gca,'XLim') get(gca,'YLim')];
	imagesc(nisy(:,:,:,frame)/255);
	axis equal, axis off, axis(current_axes)
	hold on
	%plot(pax, pay, '*r')
	hold off


	% restore original values
	nisy(ind) = nisy0(ind);

end



