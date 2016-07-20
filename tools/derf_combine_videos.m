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

	u_pat = [orig_path seq.name '/%03d.png'];                     % ground truth
%	v_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png']; % noisy
%	v_pat = [nlb3_path seq.name '_s' sigma '_pt4/deno_%03d.png']; % nlb3 pt4
%	v_pat = [nlb3_path seq.name '_s' sigma '_pt3/deno_%03d.png']; % nlb3 pt3
%	v_pat = [nlb3_path seq.name '_s' sigma '_pt2/deno_%03d.png']; % nlb3 pt2
	v_pat = [bm4d_path seq.name '_mp_s' sigma '/deno_%03d.png'];  % bm4d mp

	u = load_sequence(u_pat, seq.first, seq.last);
	v = load_sequence(v_pat, seq.first, seq.last);

%	v = combine_videos(u, v, 'vertical-half');
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

% some commands to slice video
for i = 120:-1:-120,
	slice_row = floor(size(v,1)/2) + i;
	if (size(v,3) == 3),
		slice = permute(squeeze(v(slice_row,:,:,:)),[1 3 2]);
	else
		slice = squeeze(v(slice_row,:,:));
	end

%	slice_col = floor(size(v,2)/2) + i;
%	if (size(v,3) == 3),
%		slice = permute(squeeze(v(:,slice_col,:,:)),[1 3 2]);
%	else
%		slice = squeeze(v(:,slice_col,:,:));
%	end

	imagesc(uint8(slice)',[0 255])
	%imagesc(diff(mean(slice,3),1,2)',255*[-.2 .2])
	colormap gray
	axis equal, axis off
	pause(.04)
end
