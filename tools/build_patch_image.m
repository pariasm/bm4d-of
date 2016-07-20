% A patch group is a matrix pg of size
%
%     px x px x ch x pt x np
%
% where px is the spatial  size of the patch
%       pt is the temporal size of the patch
%       ch is the number of channels
%       np is the number of patches in the group
% 
% This function transforms the patch group in to 
% an RGB image with all the patches in the group, 
% of size
%
%     pt*(px+sw)-sw x np*(px+sw)-sw x ch
%
% The temporal slices of the patch are stacked vertically,
% with separators of sw pixels. Each patch is then stacked 
% horizontally, also with separators of sw pixels.
%
% This function was written to use from patch_group.
%
function pp = build_patch_image(pp, sc, sw)

px = size(pp,1);
pt = size(pp,4);
ch = size(pp,3);
np = size(pp,5);

% separator
if nargin < 2,
	sw = 1;
	sc = 255;
end

pp = permute(pp, [2 1 4 5 3]);
pp = cat(2, pp, sc*ones(px, sw, pt, np, 3));
pp = reshape(pp, [px, pt*(px + sw), np, 3]);

pp = permute(pp, [2 1 3 4]);
pp = cat(2, pp, sc*ones(pt*(px + sw), sw, np, 3));
pp = reshape(pp, [pt*(px + sw), np*(px + sw), 3]);

pp = pp(1:pt*(px + sw) - sw, 1:np*(px + sw) - sw, :);

