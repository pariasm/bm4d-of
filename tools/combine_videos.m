% Combines two videos by spliting the frame in two
function v = combine_videos(u1, u2, mode)

% combination mode
if nargin < 3, mode = 'vertical-half'; end

if size(u1) ~= size(u2),
	error('Videos have different sizes');
end
if length(size(u1)) ~= 4,
	error('Videos must have 4 dimensions');
end

if strcmp(mode,'vertical-half') 

	half_width = floor(size(u1,2)/2);
	v = cat(2, u1(:,1:half_width,:,:),...
	           u2(:,half_width+1:end,:,:));

elseif strcmp(mode,'horizontal-half') 

	half_height = floor(size(u1,1)/2);
	v = cat(1, u1(1:half_height,:,:,:),...
	           u2(half_height+1:end,:,:,:));

elseif strcmp(mode,'horizontal-four') 

	mask = zeros(size(u1(:,:,:,1)));
	s1 = size(mask,1); s14 = floor(s1/4);
	i1 = find(mod(floor([0:s1-1]/s14),2) == 0);
	mask(i1,:,:) = 1 - mask(i1,:,:);

	mask = repmat(mask,[1,1,1,size(u1,4)]);

	v = mask.*u1 + (1-mask).*u2;

elseif strcmp(mode,'vertical-four') 

	mask = zeros(size(u1(:,:,:,1)));
	s2 = size(mask,2); s24 = floor(s2/4);
	i2 = find(mod(floor([0:s2-1]/s24),2) == 0);
	i2 = min(i2 + 30, size(mask,2));
	mask(:,i2,:) = 1 - mask(:,i2,:);
	mask = repmat(mask,[1,1,1,size(u1,4)]);

	v = mask.*u1 + (1-mask).*u2;

elseif strcmp(mode,'checkerboard-four')

	mask = zeros(size(u1(:,:,:,1)));
	s1 = size(mask,1); s14 = floor(s1/4);
	s2 = size(mask,2); s24 = floor(s2/4);
	i1 = find(mod(floor([0:s1-1]/s14),2) == 0);
	i2 = find(mod(floor([0:s2-1]/s24),2) == 0);
	mask(i1,:,:) = 1 - mask(i1,:,:);
	mask(:,i2,:) = 1 - mask(:,i2,:);

	borders = abs([diff(mask,[],1); zeros(1,size(mask,2),size(mask,3))])...
	        + abs([diff(mask,[],2)  zeros(size(mask,1),1,size(mask,3))]);

	mask    = repmat(mask   ,[1,1,1,size(u1,4)]);
	borders = repmat(borders,[1,1,1,size(u1,4)]);

	v = mask.*u1 + (1-mask).*u2;

	ib = find(borders ~= 0);
	v(ib) = 0;

else
	error('Unknown combination mode %s', mode);
end
