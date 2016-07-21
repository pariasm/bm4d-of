root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d-dct/results/';
root_folder =  '/home/pariasm/Work/denoising/projects/video_nlbayes3d-dct/results/';

stage1_folder_pos = [root_folder 'table_1np_vs_2np_img-color_1ctrY_1thN_1px4_2px4/'];
stage1_folder_neg = [root_folder 'table_1np_vs_2np_img-color_1ctrY_1thN_1px8_2px4/'];
%stage1_folder_neg = [root_folder 'table_1np_vs_2np_img-color_1thN_2thN_1px8_2px8/'];

method1 = 'WIn4WIE4';
method2 = 'WIn8WIE4';

% stage 1 ---------------------------------------------------
nsims1_p = [10 40 100];
nsims1_n = [10 40 100];
nsims2_p = [6 10 20 40 60 100];
nsims2_n = [6 10 20 40 60 100];
seqs = {...
'alley',...
'building1',...
'computer',...
'dice',...
'girl',...
'street1',...
'traffic',...
'trees'};

sigmas = [10,40];

eff_p = zeros(length(nsims1_p),length(nsims2_p),length(sigmas));
ebb_p = zeros(length(nsims1_p),length(nsims2_p),length(sigmas));

eff_n = zeros(length(nsims1_n),length(nsims2_n),length(sigmas));
ebb_n = zeros(length(nsims1_n),length(nsims2_n),length(sigmas));


% indices to be removed
idx_n1_n = [];
idx_n1_p = [];

idx_n2_n = [];
idx_n2_p = [];

fig_i = 10;
for sidx = 1:length(sigmas), s = sigmas(sidx);

	% load final neg
	ff_n = nan(length(nsims1_n),length(nsims2_n),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_neg 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp(:,idx_n2_n) = [];
		tmp(idx_n1_n,:) = [];
		ff_n(:,:,i) = tmp;
	end

	% load basic neg
	bb_n = nan(length(nsims1_n),length(nsims2_n),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_neg 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp(:,idx_n2_n) = [];
		tmp(idx_n1_n,:) = [];
		bb_n(:,:,i) = tmp;
	end

	% load final pos
	ff_p = nan(length(nsims1_p),length(nsims2_p),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_pos 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		% remove indices that have not been computed
		tmp(:,idx_n2_p) = [];
		tmp(idx_n1_p,:) = [];
		ff_p(:,:,i) = tmp;
	end

	% load basic pos
	bb_p = nan(length(nsims1_p),length(nsims2_p),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_pos 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp(:,idx_n2_p) = [];
		tmp(idx_n1_p,:) = [];
		bb_p(:,:,i) = tmp;
	end

	% ------------------
	% average final psnr
	eff_p(:,:,sidx) = mean(ff_p,3);
	eff_n(:,:,sidx) = mean(ff_n,3);

	figure(fig_i), clf,
	fig_i = fig_i + 1;

	m_p = max(max(eff_p(:,:,sidx)));
	m_n = max(max(eff_n(:,:,sidx)));
	m = max(m_n, m_p);

	hold on
	color_step = 1/length(nsims1_p);
	for ii = length(nsims1_p):-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(nsims2_p, eff_p(ii,:,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','o','MarkerSize',5,'MarkerFaceColor',color)
		axis([0 100 m-1 m+.2])
		grid on
	end
	color_step = 1/length(nsims1_n);
	for ii = length(nsims1_n):-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(nsims2_n, eff_n(ii,:,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','s','MarkerSize',5,'MarkerFaceColor',color)
		axis([0 100 m-0.8 m+.2])
		grid on
	end
	hold off
	h_leg = legend(...
	               '$P1 - n_1 = 100$',...
	               '$P1 - n_1 =  40$',...
	               '$P1 - n_1 =  10$',...
	               '$P2 - n_1 = 100$',...
	               '$P2 - n_1 =  40$',...
	               '$P2 - n_1 =  10$',...
	               'Location','SouthEast');
	set(h_leg,'Interpreter','latex');
	h_xlab = xlabel('$n_2$','Interpreter','latex');
	h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');

	% NewCenturySchlbk, AvantGarde, Helvetica
	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	box on

%	% leave only integer labels, so that final PsNR and basic PSNR figures have
%	% same width
%	ylls = get(gca,'yticklabel');
%	for ll = 1:size(ylls,1),
%		if mod(str2num(ylls(ll,:)),1) ~= 0,
%			ylls(ll,:) = ' ' * ones(1,size(ylls,2));
%		end
%	end
%	set(gca,'yticklabel', ylls(:,1:2));

%	print(gcf, '-depsc2', sprintf('%s-vs-%s_fpsnr_1np-np2-curves_s%02d_average',method1, method2, s));
%	print(gcf, '-dpng',   sprintf('%s-vs-%s_fpsnr_1np-np2-curves_s%02d_average',method1, method2, s));

	% ------------------
	% average basic psnr
	ebb_p(:,:,sidx) = mean(bb_p,3);
	ebb_n(:,:,sidx) = mean(bb_n,3);


end



