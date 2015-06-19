% Computes the evolution of the frame psnr and mser.
function [psnr, mser] = psnr_curves(v, u)

if size(u) ~= size(v),
	error('Videos have different sizes');
end
if length(size(u)) ~= 4,
	error('Videos must have 4 dimensions');
end

imax = max(v(:));
if imax > 40, imax = 255; else, imax = 1; end

frames = size(v,4);

psnr = zeros(frames,1);
mser = zeros(frames,1);

for i = 1:frames,
	u_i = u(:,:,:,i);
	v_i = v(:,:,:,i);
	mser(i) = mean((u_i(:) - v_i(:)).^2);
	psnr(i) = 10*log10(imax^2/mser(i));
end

