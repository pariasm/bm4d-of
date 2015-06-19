% combines different denoising results for in the same video,
% for visualization purposes

% path to data
base_path = '/media/pariasm/tera/funes/denoising/';
bm4d_path = [base_path 'algos/results/VBM4D/'];
orig_path = [base_path 'data/middlebury/'];
nlb3_path = [base_path 'projects/video_nlbayes3d/results/vnlbayes/table_1ps5_2ps5_2wx37/'];
nlb2_path = [base_path 'projects/video_nlbayes/results/vnlbayes/table_1ps5_2ps5_2wx37_2np93/'];

% sequence data
seqnames = {'Army','DogDance','Evergreen','Mequon','Walking'};
first    = {0     , 0        , 0         , 0      , 0       };
last     = {7     , 7        , 7         , 7      , 7       };

sequences = struct('name',seqnames,'first',first,'last',last);
clear seqnames first last

% parameters
sigma = '40';

addpath '../../../../time_filter/src/';

for seqi = [1],%1:length(sequences),

	seq = sequences(seqi);

	v_pat = [orig_path seq.name '/i%04d.png'];                     % ground truth
%	v_pat = [nlb3_path seq.name '_s' sigma '_pt4/nisy_%03d.png']; % noisy
	u_pat = [nlb3_path seq.name '_s' sigma '_pt4/deno_%03d.png']; % nlb3 pt4
%	u_pat = [nlb3_path seq.name '_s' sigma '_pt3/deno_%03d.png']; % nlb3 pt3
%  u_pat = [nlb3_path seq.name '_s' sigma '_pt2/deno_%03d.png']; % nlb3 pt2
%	u_pat = [nlb2_path seq.name '_s' sigma '_tr2/deno_%03d.png']; % nlb3 pt2
%	v_pat = [nlb2_path seq.name '_s' sigma '_tr2/deno_%03d.png']; % nlb3 pt2
%	v_pat = [bm4d_path seq.name '_s' sigma '/deno_%03d.png'];  % bm4d mp

	u = load_sequence(u_pat, seq.first, seq.last);
	v = load_sequence(v_pat, seq.first, seq.last);

	v = combine_videos(u, v, 'checkerboard-four');
%	save_sequence(v, '/tmp/combi_%03d.png', seq.first);
	o_dir = ['combined_videos/' seq.name '_s' sigma];
	mkdir(o_dir);
	save_sequence(v, [o_dir '/vnlb_pt4-orig_%03d.png'], seq.first);


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

% % some commands to slice video
% for i = 120:-1:-120,
% 	slice_row = floor(size(v,1)/2) + i;
% 	if (size(v,3) == 3),
% 		slice = permute(squeeze(v(slice_row,:,:,:)),[1 3 2]);
% 	else
% 		slice = squeeze(v(slice_row,:,:));
% 	end
% 
% %	slice_col = floor(size(v,2)/2) + i;
% %	if (size(v,3) == 3),
% %		slice = permute(squeeze(v(:,slice_col,:,:)),[1 3 2]);
% %	else
% %		slice = squeeze(v(:,slice_col,:,:));
% %	end
% 
% 	%imagesc(uint8(slice),[0 255])
% 	imagesc(diff(mean(slice,3),1,2),255*[-.2 .2])
% 	colormap gray
% 	axis equal, axis off
% 	pause(.04)
% end
