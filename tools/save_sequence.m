function save_sequence(v, inpat, from)

if nargin < 3,
	from = 1;
end

% load first frame to get dimensions, channels
sz = [size(v,1) size(v,2)];
ch =  size(v,3);
nf =  size(v,4);

% load sequence 
disp(sprintf('Writing %s', inpat))
for i = 1:nf,

	filename = sprintf(inpat,from + i - 1);
	imwrite(uint8(squeeze(v(:,:,:,i))),filename);

end
