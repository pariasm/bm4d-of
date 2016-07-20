% plot the mser and psnr curves for some 
% sequences of the derf dataset

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37_2r16_2np160/'];
bm4d_path = [base_path 'algos/results/VBM4D/'];
orig_path = [base_path 'data/derf/'];

% sequence data
seqnames = {'bus','coastguard','foreman','tennis'};
first    = {1    , 1          , 1       , 1      };
last     = {150  , 300        , 300     , 150    };

sequences = struct('name',seqnames,'first',first,'last',last);
clear seqnames first last

% parameters
sigma = '40';
%pt = '4';

addpath '../../../../time_filter/src/';

for seqi = 1:length(sequences),

	seq = sequences(seqi);

	% load ground truth
	u = load_sequence([orig_path seq.name '/%03d.png'], seq.first, seq.last);

	% load noisy sequence
	seq_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png'];
	v = load_sequence(seq_pat, seq.first, seq.last);
	[psnrn, msern] = psnr_curves(u, v);

	% load denoised with pt4
	seq_pat = [nlb3_path seq.name '_s' sigma '_pt4/deno_%03d.png'];
	v = load_sequence(seq_pat, seq.first, seq.last);
	[psnr4, mser4] = psnr_curves(u, v);

	% load denoised with pt3
	seq_pat = [nlb3_path seq.name '_s' sigma '_pt3/deno_%03d.png'];
	v = load_sequence(seq_pat, seq.first, seq.last);
	[psnr3, mser3] = psnr_curves(u, v);

	% load denoised with pt2
	seq_pat = [nlb3_path seq.name '_s' sigma '_pt2/deno_%03d.png'];
	v = load_sequence(seq_pat, seq.first, seq.last);
	[psnr2, mser2] = psnr_curves(u, v);

	% load result of v-bm4d
	%seq_pat = [bm4d_path seq.name '_mp_s' sigma '/deno_%03d.png'];
	%v = load_sequence(seq_pat, seq.first, seq.last);
	[psnr_bm4d, mser_bm4d] = psnr_curves(u, v);


	% plot results
	figure(1),
	ff = [seq.first:seq.last];
	plot(ff, sqrt(msern    ) , ...
	     ff, sqrt(mser4    ) , ...
	     ff, sqrt(mser3    ) , ...
	     ff, sqrt(mser2    ) , ...
	     ff, sqrt(mser_bm4d) )
	legend('noisy','pt4','pt3','pt2','V-BM4D-mp')
	title(seq.name)

	figure(2),
	ff = [seq.first:seq.last];
	plot(ff, psnrn     , ...
	     ff, psnr4     , ...
	     ff, psnr3     , ...
	     ff, psnr2     , ...
	     ff, psnr_bm4d )
	legend('noisy','pt4','pt3','pt2','V-BM4D-mp')
	title(seq.name)


	keyboard

end

