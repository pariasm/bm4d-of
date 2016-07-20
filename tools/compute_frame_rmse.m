function rmse = compute_frame_rmse(vin, gt)

h  = size(vin,1);
w  = size(vin,2);
nc = size(vin,3);
nf = size(vin,4);

rmse = zeros(nf,1);

for f = 1:nf,

	tmp = vin(:,:,:,f) - gt(:,:,:,f);
	rmse(f) = sqrt(mean(tmp(:).^2));

end


