% downscale by a factor of 2 using lanczos3 kernel
function down = lanczos3_down(im);

% -----x---x---x---x---x---x---x---x---x---x----- 2n high res signal
% -------o-------o-------o-------o-------o------- n low res

% input size
h = size(im,1);
w = size(im,2);

k = 0.5 * lanczos3_kernel(0.5* (0.5 + [-6:5])); k = k/sum(k(:));
p = [repmat(im(:,1),1,5), im , repmat(im (:,w),1,6)]; tmp = conv2(p, k, 'valid');
down1 = tmp(:,1:2:end);

p = [repmat(down1(1,:),5,1); down1 ; repmat(down1(h,:),6,1)]; tmp = conv2(p, k', 'valid');
down = tmp(1:2:end,:);

