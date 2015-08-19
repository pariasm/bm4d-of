% combines different denoising results for in the same video,
% for visualization purposes

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
bm4d_path = [base_path 'algos/results/VBM4D/'];
orig_path = [base_path 'data/derf/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37_2r16_2np160/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps7_2ps5_2wx45_r16/'];
nlb3_colo_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_tip/'];
nlb3_mono_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps7_2ps7_tip/'];
nisy_path = bm4d_path;


% generated tiles ------------------- %
% tennis       050,090,030 150,150,30 %
% tennis       040,150,070 150,150,30 %
% bus          060,020,020 150,150,30 %
% bus          110,150,010 150,150,30 %
% stefan_mono  090,180,265 150,150,30 %
% stefan_mono  060,130,005 150,150,30 %
% mobile_mono  060,060,080 100,100,30 %
% mobile_mono  188,188,065 100,100,30 %
% ----------------------------------- %

% sequence data
sequences = struct([]);

sequences(1).name = 'football_mono';
sequences(1).first = 1;
sequences(1).last = 260;
sequences(1).channels = 1;

sequences(2).name = 'mobile_mono';
sequences(2).first = 1;
sequences(2).last = 300;
sequences(2).channels = 1;

sequences(3).name = 'stefan_mono';
sequences(3).first = 1;
sequences(3).last = 300;
sequences(3).channels = 1;

sequences(4).name = 'tennis_mono';
sequences(4).first = 1;
sequences(4).last = 150;
sequences(4).channels = 1;

sequences(5).name = 'bus';
sequences(5).first = 1;
sequences(5).last = 150;
sequences(5).channels = 3;

sequences(6).name = 'coastguard';
sequences(6).first = 1;
sequences(6).last = 300;
sequences(6).channels = 3;

sequences(7).name = 'foreman';
sequences(7).first = 1;
sequences(7).last = 300;
sequences(7).channels = 3;

sequences(8).name = 'tennis';
sequences(8).first = 1;
sequences(8).last = 150;
sequences(8).channels = 3;


tiles = struct([]);

tiles(1).seq = sequences(8);
tiles(1).tile_b = [050,090,030];
tiles(1).tile_e = tiles(1).tile_b + [150,150, 30];

tiles(2).seq = sequences(8);
tiles(2).tile_b = [040,150,070];
tiles(2).tile_e = tiles(2).tile_b + [150,150, 30];

tiles(3).seq = sequences(5);
tiles(3).tile_b = [060,020,020];
tiles(3).tile_e = tiles(3).tile_b + [150,150, 30];

tiles(4).seq = sequences(5);
tiles(4).tile_b = [110,150,010];
tiles(4).tile_e = tiles(4).tile_b + [150,150, 30];

tiles(5).seq = sequences(3);
tiles(5).tile_b = [090,180,265];
tiles(5).tile_e = tiles(5).tile_b + [150,150, 30];

tiles(6).seq = sequences(3);
tiles(6).tile_b = [060,130,005];
tiles(6).tile_e = tiles(6).tile_b + [150,150, 30];

tiles(7).seq = sequences(2);
tiles(7).tile_b = [060,060,080];
tiles(7).tile_e = tiles(7).tile_b + [100,100, 30];

tiles(8).seq = sequences(2);
tiles(8).tile_b = [188,188,065];
tiles(8).tile_e = tiles(8).tile_b + [100,100, 30];


%addpath '../../time_filter/src/';


sigmas = {'10', '20', '40'};

%for vnlb = ['O','N','B','4'],
for vnlb = ['B','4'],
for sigmai = 1:3, sigma = sigmas{sigmai};
for t = 1:length(tiles),

	seq = tiles(t).seq;

	switch (vnlb)
	case 'B',
		if seq.channels == 1,
			u_pat = [ nlb3_mono_path seq.name '_s' sigma '_pt4/deno_%03d.png']; % nlb3 pt4
		else
			u_pat = [ nlb3_colo_path seq.name '_s' sigma '_pt4/deno_%03d.png']; % nlb3 pt4
		end
	case '4', u_pat = [bm4d_path seq.name '_mp_s' sigma '/deno_%03d.png'];  % bm4d mp
	case 'O', u_pat = [orig_path seq.name '/%03d.png'];                     % ground truth
	case 'N', u_pat = [nisy_path seq.name '_mp_s' sigma '/nisy_%03d.png'];  % noisy
	end

	u = load_sequence(u_pat, seq.first, seq.last);

	% split video in parts
	sz = size(u);
	tile_sz = [floor(sz(1)/3) floor(sz(2)/4) sz(3) sz(4)];

%	for i = 0:floor(sz(1)/tile_sz(1))-2,
%	for j = 0:floor(sz(2)/tile_sz(2))-2,
%		save_sequence(u(i*tile_sz(1) + floor(tile_sz(1)/2) + [1:tile_sz(1)], ...
%		                j*tile_sz(2) + floor(tile_sz(2)/2) + [1:tile_sz(2)], :, :), ...
%		              sprintf('/tmp/tile%d_%d_%%03d.png', i, j));
%	end
%	end

	switch (vnlb)
%	case 'B', tile_pat = [nlb3_path seq.name '_s' sigma '_pt4/tile%03d-%03d-%03d_%%03d.png']; % nlb3 pt4
%	case '4', tile_pat = [bm4d_path seq.name '_mp_s' sigma  '/tile%03d-%03d-%03d_%%03d.png']; % bm4d mp
%	case 'O', tile_pat = [orig_path seq.name '/tile%03d-%03d-%03d_%%03d.png'];                     % ground truth
%	case 'N', tile_pat = [nisy_path seq.name '_mp_s' sigma '/nisy_tile%03d-%03d-%03d_%%03d.png']; % noisy
	case 'B', tile_pat = [seq.name '_s' sigma '_pt4_tile%03d-%03d-%03d_%%03d.png']; % nlb3 pt4
	case '4', tile_pat = [seq.name '_mp_s' sigma  '_tile%03d-%03d-%03d_%%03d.png']; % bm4d mp
	case 'O', tile_pat = [seq.name '_tile%03d-%03d-%03d_%%03d.png'];                     % ground truth
	case 'N', tile_pat = [seq.name '_mp_s' sigma '_nisy_tile%03d-%03d-%03d_%%03d.png']; % noisy
	end

	tile_i = [tiles(t).tile_b(1) : tiles(t).tile_e(1)];
	tile_j = [tiles(t).tile_b(2) : tiles(t).tile_e(2)];
	tile_t = [tiles(t).tile_b(3) : tiles(t).tile_e(3)];
	save_sequence(u(tile_i, tile_j, :, tile_t), ...
	              sprintf(tile_pat, tile_i(1), tile_j(1), tile_t(1)));


	%save_sequence(v, '/tmp/combi_%03d.png', seq.first);


%	% plot results
%	figure(1),
%	ff = [seq.first:seq.last];
%	plot(ff, sqrt(msern    ) , ...
%	     ff, sqrt(mser4    ) , ...
%	     ff, sqrt(mser3    ) , ...
%	     ff, sqrt(mser2    ) , ...
%	     ff, sqrt(mser_bm4d) )
%	legend('noisy','pt4','pt3','pt2','V-BM4D-mp')
%	title(seq.name)
%
%	figure(2),
%	ff = [seq.first:seq.last];
%	plot(ff, psnrn     , ...
%	     ff, psnr4     , ...
%	     ff, psnr3     , ...
%	     ff, psnr2     , ...
%	     ff, psnr_bm4d )
%	legend('noisy','pt4','pt3','pt2','V-BM4D-mp')
%	title(seq.name)

end
end
end
