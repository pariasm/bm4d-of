% combines different denoising results for in the same video,
% for visualization purposes

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
bm4d_path = [base_path 'algos/results/VBM4D/'];
orig_path = [base_path 'data/derf/'];
%nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37_2r16_2np160/'];
nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps7_2ps5_2wx45_r16/'];

% sequence data
seqnames = {'football_mono','mobile_mono','stefan_mono','tennis_mono'};
first    = {1              , 1           , 1           , 1           };
last     = {260            , 300         , 300         , 150         };
%seqnames = {'bus','coastguard','foreman','tennis'};
%first    = {1    , 1          , 1       , 1      };
%last     = {150  , 300        , 300     , 150    };

sequences = struct('name',seqnames,'first',first,'last',last);
clear seqnames first last

% parameters
sigma = '10';

addpath '../../time_filter/src/';

for seqi = [2],%1:length(sequences),

	seq = sequences(seqi);

	for pt = 2:7,
		if pt < 5,
			u_pat = [nlb3_path seq.name '_s' sigma '_pt' num2str(pt) '/deno_%03d.png']; % nlb3 pt4
			out_name = ['_vnlb_pt' num2str(pt) '_s' sigma];
		elseif pt == 5,
			u_pat = [bm4d_path seq.name '_mp_s' sigma '/deno_%03d.png'];  % bm4d mp
			out_name = ['_bm4d' '_s' sigma];
		elseif pt == 6,
			u_pat = [orig_path seq.name '/%03d.png'];
			out_name = '_orig';
		elseif pt == 7,
			u_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png']; % noisy
			out_name = ['_nisy' '_s' sigma];
		end

		u = load_sequence(u_pat, seq.first, seq.last);

		% slice
		slice = squeeze(u(220,40:180,:,80:220))';

		% linear scaling
		slice = uint8(255*(slice - 80)/80);

		% write
		imwrite(slice, ['slice_' seq.name out_name  '_row220_col040-180_fra080-220.png']);

	end
end

% now save the error
for seqi = [2],%1:length(sequences),

	seq = sequences(seqi);

	% load original sequence
	u_pat = [orig_path seq.name '/%03d.png'];
	u = load_sequence(u_pat, seq.first, seq.last);
	slice_orig = squeeze(u(220,40:180,:,80:220))';

	for pt = 2:5,
		if pt < 5,
			u_pat = [nlb3_path seq.name '_s' sigma '_pt' num2str(pt) '/deno_%03d.png']; % nlb3 pt4
			out_name = ['_diff_vnlb_pt' num2str(pt) '_s' sigma];
		else pt == 5,
			u_pat = [bm4d_path seq.name '_mp_s' sigma '/deno_%03d.png'];  % bm4d mp
			out_name = ['_diff_bm4d' '_s' sigma];
		end

		u = load_sequence(u_pat, seq.first, seq.last);

		% slice
		slice = squeeze(u(220,40:180,:,80:220))';

		% linear scaling
		slice = uint8(255*(((slice - slice_orig + 120) - 80)/80));

		% write
		imwrite(slice, ['slice_' seq.name out_name  '_row220_col040-180_fra080-220.png']);

	end
end
