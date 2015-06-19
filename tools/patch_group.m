% loads a video and allows the user to select a pixel
% for the selected pixel, it displays the locations 
% of the most similar patches

clear all


% parameters
prms.wx = '37';
prms.wt = '2';
prms.px = '19';
prms.pt = '4';
prms.np = '120';
sigma = '40';

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
nsy_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png']; % noisy
bsc_pat = [nlb3_path seq.name '_s' sigma '_pt4/bsic_%03d.png']; % nlb3 pt2

command = [pgbin ...
           ' -i ' ori_pat ...
           ' -b ' ori_pat ...
           ' -f ' num2str(seq.first) ...
           ' -l ' num2str(seq.last) ...
           ' -sigma ' sigma ...
			  ' -px ' prms.px ...
			  ' -pt ' prms.pt ...
			  ' -wx ' prms.wx ...
			  ' -wt ' prms.wt ...
			  ' -np ' prms.np   ];

% -PAx 100 -PAy 200 -PAt 3  > patches

% load
u = load_sequence(ori_pat, seq.first, seq.last);
u0 = u;

% show sequence
goon = 1;
frame = 1;
pax = -1;
pay = -1;
pat = -1;

imagesc(u(:,:,:,frame)/255);
axis equal, axis off

while (goon)

	key = input('','s');

	switch (key),
	case 'f',
		frame = min(frame + 1, size(u,4));
	case 'b',
		frame = max(frame - 1, 1);
	case 's',
		[pax,pay] = ginput(1)
		pax = round(pax);
		pay = round(pay);
		pat = frame;

		% compute patch group
		cmd = [command ...
				' -PAx ' num2str(pax - 1) ... 
				' -PAy ' num2str(pay - 1) ... 
				' -PAt ' num2str(frame - 1) ... 
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

		ind = sub2ind(size(u),patches.coords(:,2), ...
		                      patches.coords(:,1), ...
		                      ones(size(patches.coords,1),1), ...
									 patches.coords(:,3));
		u(ind) = 255;

		patches.vals = zeros(19,19,3,4,length(patches.distas));
		for i = 1:length(patches.distas),
			patches.vals(:,:,:,:,i) = u(patches.coords(i,2):patches.coords(i,2)+19-1,...
			                            patches.coords(i,1):patches.coords(i,1)+19-1,:,...
			                            patches.coords(i,3):patches.coords(i,3)+ 4-1);
		end


	case 'q'
		goon = 0;
	otherwise,
	end

	imagesc(u(:,:,:,frame)/255);
	axis equal, axis off
	hold on
	plot(pax, pay, '*r')
	hold off

end



