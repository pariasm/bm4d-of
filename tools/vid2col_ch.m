% Video version of im2col_ch. Rearranges video blocks into columns.
% For a multichannel video, the vectorized blocks of each channel are staked
% vertically in the output.
%
% USAGE: patches = vid2col_ch(vid, [ph pw pf], mode)
%
%  -> vid      : input video (height x width x channels x frames)
%  -> ph,pw,pf : patch size
%  -> mode     : either 'distinct' of 'sliding'
%
%  <- patches : output patches (ph pw ch pf x n)
function patches = vid2col_ch(vid, psz, mode)

	if nargin < 3,
		mode = 'sliding';
	end

	ph = psz(1);
	pw = psz(2);
	pf = psz(3);
	[h w c f] = size(vid);

	if pf <= f,

		% patch dimension (single channel)
		pdim = prod(psz);
		pdim2D = prod(psz(1:2));
		pdim2Dc = pdim2D*c;

		% get the number of patches that fit in each dimension
		if strcat(mode,'sliding')
			time_step = 1;
			nt = f - pf + 1;
			nx = w - pw + 1;
			ny = h - ph + 1;
		else
			time_step = pf;
			nt = ceil(f/pf);
			nx = ceil(w/pw);
			ny = ceil(h/ph);
		end
		nxy = nx*ny;

		% get the number of patches per frame from first frame, first channel
		patches = zeros(c*pdim, nx*ny*nt);


		if strcat(mode,'sliding')

			for t = 1:f,

				tmin = max(t - pf + 1,1);
				tmax = min(f - pf + 1,t);

				for k = 0:c-1,
					patches2D = im2col(vid(:,:,k+1,t),[ph pw], mode);

					for t1 = tmax:-1:tmin,
						patches((t-t1)*pdim2Dc + k*pdim2D + [1:pdim2D], (t1-1)*nxy + [1:nxy]) = patches2D;
					end
				end
			end

		else

			error('This mode is not yet implemented')

		end

	else
		patches = [];
	end

end



