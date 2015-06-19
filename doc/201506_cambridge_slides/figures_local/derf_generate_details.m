% combines different denoising results for in the same video,
% for visualization purposes

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
bm4d_path = [base_path 'algos/results/VBM4D/'];
orig_path = [base_path 'data/derf/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37_2r16_2np160/'];
nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps7_2ps5_2wx45_r16/'];


% generated tiles ------------------- %
% tennis       050-090-030 150+150+30 %
% tennis       040-150-070 150+150+30 %
% bus          060-020-020 150+150+30 %
% bus          110-150-010 150+150+30 %
% stefan_mono  090-180-265 150+150+30 %
% stefan_mono  060-130-005 150+150+30 %
% mobile_mono  060-060-080 100+100+30 %
% mobile_mono  188-188-065 100+100+30 %
% ----------------------------------- %

ti = 060;
tj = 060;
tt = 080;
tile_i = [ti:ti + 100];
tile_j = [tj:tj + 100];
tile_t = [tt:tt +  30];

addpath '../../../tools/';

% sequence data
seqnames = {'football_mono','mobile_mono','stefan_mono','tennis_mono'};
first    = {1              , 1           , 1           , 1           };
last     = {260            , 300         , 300         , 150         };
%seqnames = {'bus','coastguard','foreman','tennis'};
%first    = {1    , 1          , 1       , 1      };
%last     = {150  , 300        , 300     , 150    };

sequences = struct('name',seqnames,'first',first,'last',last);
clear seqnames first last

sigmas = {'10', '20', '40'};

%for vnlb = ['O','N','B','4'],
for vnlb = ['N'],
for sigmai = 1:3, sigma = sigmas{sigmai};
for seqi = [2],

	seq = sequences(seqi);

	switch (vnlb)
	case 'B', u_pat = [nlb3_path seq.name '_s' sigma '_pt4/deno_%03d.png']; % nlb3 pt4
	case '4', u_pat = [bm4d_path seq.name '_mp_s' sigma '/deno_%03d.png'];  % bm4d mp
	case 'O', u_pat = [orig_path seq.name '/%03d.png'];                     % ground truth
	case 'N', u_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png']; % noisy
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
	case 'B', tile_pat = [nlb3_path seq.name '_s' sigma '_pt4/tile%03d-%03d-%03d_%%03d.png']; % nlb3 pt4
	case '4', tile_pat = [bm4d_path seq.name '_mp_s' sigma  '/tile%03d-%03d-%03d_%%03d.png']; % bm4d mp
	case 'O', tile_pat = [orig_path seq.name '/tile%03d-%03d-%03d_%%03d.png'];                     % ground truth
	case 'N', tile_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_tile%03d-%03d-%03d_%%03d.png']; % noisy
	end

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

	break

end
end
end
