% combines different denoising results for in the same video,
% for visualization purposes

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
bm4d_path = [base_path 'algos/results/VBM4D/'];
orig_path = [base_path 'data/derf/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37_2r16_2np160/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps7_2ps5_2wx45_r16/'];
nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/tip_table_1ps7_2ps7-neg/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/tip_table_1ps5_2ps5-neg/'];


% generated tiles ------------------- %

%seqnames = {'football_mono','mobile_mono','stefan_mono','tennis_mono'};
%first    = {1              , 1           , 1           , 1           };
%last     = {260            , 300         , 300         , 150         };
%seqnames = {'bus','coastguard','foreman','tennis'};
%first    = {1    , 1          , 1       , 1      };
%last     = {150  , 300        , 300     , 150    };

% color

%---% seq.name = 'tennis'; seq.first = 1; seq.last = 150;
%---% ti=050; tj=090; tt=030; wi=150; wj=150; wt=30;
%---% 
%---% seq.name = 'tennis'; seq.first = 1; seq.last = 150;
%---% ti=040; tj=150; tt=070; wi=150; wj=150; wt=30;
%---% 
%---% seq.name = 'bus'; seq.first = 1; seq.last = 150;
%---% ti=060; tj=020; tt=020; wi=150; wj=150; wt=30;
%---% 
%---% seq.name = 'bus'; seq.first = 1; seq.last = 150;
%---% ti=110; tj=150; tt=010; wi=150; wj=150; wt=30;

% grayscale

%---% seq.name = 'stefan_mono'; seq.first = 1; seq.last = 300;
%---% ti=090; tj=180; tt=265; wi=150; wj=150; wt=30;
%---% 
%---% seq.name = 'stefan_mono'; seq.first = 1; seq.last = 300;
%---% ti=060; tj=130; tt=005; wi=150; wj=150; wt=30;
%---% 
%---% seq.name = 'mobile_mono'; seq.first = 1; seq.last = 300;
%---% ti=060; tj=060; tt=080; wi=100; wj=100; wt=30;
%---% 
seq.name = 'mobile_mono'; seq.first = 1; seq.last = 300;
ti=188; tj=188; tt=065; wi=100; wj=100; wt=30;
% ----------------------------------- %


tile_i = [ti:ti + wi];
tile_j = [tj:tj + wj];
tile_t = [tt:tt + wt];

addpath '../../../tools/';


sigmas = {'10', '20', '40'};

%for vnlb = ['O','N','B','4'],
for vnlb = ['B'],
for sigmai = 1:3, sigma = sigmas{sigmai};

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
