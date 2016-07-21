root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d-dct/results/';
root_folder =  '/home/pariasm/Work/denoising/projects/video_nlbayes3d-dct/results/';

%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1lzd_2wie_1px8_2px8/'];
stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1ctrY_1thN_1px8_2px8/'];
stage1_folder_m2 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1ctrY_1thY_1px8_2px8/'];
%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1thY_1mu3_2muN_1px8_2px8/'];
%stage1_folder_m2 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1ctrY_1thY_b2f_1px8_2px8/'];
%stage1_folder_m2 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1thN_1mu3_2muN_1px8_2px8/'];
%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1softB_1px8_2px8/'];
%stage1_folder_m2 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1softB_b2f_1px8_2px8/'];
%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1linrT_1px8_2px8/'];
%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1hardT_1px8_2px8/'];
%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1softB_1px8_2px8/'];
%stage1_folder_m2 = [root_folder 'table_1np_vs_2np_vs_1b_img-color_1hardT_1px8_2px8/'];
%stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_2b_img-color_1thB_2thB_1px8_2px8/'];
%stage1_folder_m2 = [root_folder 'table_1np_vs_2np_vs_2b_img-color_1wiN_2wie_1px8_2px8/'];

method1 = 'WIE1neg';
method2 = 'WIE1pos';

nsims1_m1 = [10 20 30 40 100];
nsims1_m2 = [10 20 30 40 100];
nsims2_m1 = [6 10 20 40 60];
nsims2_m2 = [6 10 20 40 60];
%b1s_m1 = [0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4];
%b1s_m2 = [0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4];
%b1s_m1 = [0.70 : 0.05 : 1.45];
b1s_m1 = [0.70 : 0.05 : 1.45];
b1s_m2 = [0.70 : 0.05 : 1.45];
%b1s_m2 = [0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4];
%b1s_m1 = [0.4 0.8 1.0 1.2 1.6 2.0];
%b1s_m2 = [0.4 0.8 1.0 1.2 1.6 2.0];
%b1s_m1 = [0.4 : 0.2 : 2.2];
%b1s_m2 = [0.4 : 0.2 : 2.2];
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

eff_m1 = zeros(length(nsims2_m1),length(b1s_m1),length(nsims1_m1),length(sigmas));
ebb_m1 = zeros(length(nsims2_m1),length(b1s_m1),length(nsims1_m1),length(sigmas));
                                                              
eff_m2 = zeros(length(nsims2_m2),length(b1s_m2),length(nsims1_m2),length(sigmas));
ebb_m2 = zeros(length(nsims2_m2),length(b1s_m2),length(nsims1_m2),length(sigmas));


% indices to be removed
idx_n1_m2 = [];
idx_n1_m1 = [];

idx_n2_m2 = [];
idx_n2_m1 = [];

fig_i = 10;
for sidx = 1:length(sigmas), s = sigmas(sidx);

	% load final m2
	ff_m2 = nan(length(nsims2_m2),length(b1s_m2),length(nsims1_m2),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_m2 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp(:,idx_n2_m2) = [];
		tmp(idx_n1_m2,:) = [];
		for j = 1:length(nsims1_m2),
			ff_m2(:,:,j,i) = tmp([1:length(nsims2_m2)] + length(nsims2_m2)*(j-1),:);
		end
	end

	% load basic m2
	bb_m2 = nan(length(nsims2_m2),length(b1s_m2),length(nsims1_m2),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_m2 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp(:,idx_n2_m2) = [];
		tmp(idx_n1_m2,:) = [];
		for j = 1:length(nsims1_m2),
			bb_m2(:,:,j,i) = tmp([1:length(nsims2_m2)] + length(nsims2_m2)*(j-1),:);
		end
	end

	% load final m1
	ff_m1 = nan(length(nsims2_m1),length(b1s_m1),length(nsims1_m1),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_m1 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		% remove indices that have not been computed
		tmp(:,idx_n2_m1) = [];
		tmp(idx_n1_m1,:) = [];
		for j = 1:length(nsims1_m1),
			ff_m1(:,:,j,i) = tmp([1:length(nsims2_m1)] + length(nsims2_m1)*(j-1),:);
		end
	end

	% load basic m1
	bb_m1 = nan(length(nsims2_m1),length(b1s_m1),length(nsims1_m1),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_m1 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp(:,idx_n2_m1) = [];
		tmp(idx_n1_m1,:) = [];
		for j = 1:length(nsims1_m1),
			bb_m1(:,:,j,i) = tmp([1:length(nsims2_m1)] + length(nsims2_m1)*(j-1),:);
		end
	end

	if 1,
	% ------------------
	% average final psnr
	eff_m1(:,:,:,sidx) = mean(ff_m1,4);
	eff_m2(:,:,:,sidx) = mean(ff_m2,4);


	m_m1 = max(max(max(eff_m1(:,:,:,sidx))));
	m_m2 = max(max(max(eff_m2(:,:,:,sidx))));
	m = max(m_m2, m_m1);

	for jj = 1:length(nsims1_m1),

		if nsims1_m1(jj) == 20 | nsims1_m1(jj) == 30,
			continue;
		end


		figure(fig_i), clf,
		fig_i = fig_i + 1;

		hold on
		color_step = 1/length(nsims2_m1);
		for ii = length(nsims2_m1):-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(b1s_m1, eff_m1(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2)
		end
		for ii = length(nsims2_m1):-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(b1s_m1, eff_m1(ii,:,jj,sidx),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
		end
		color_step = 1/length(nsims2_m2);
		for ii = length(nsims2_m2):-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(b1s_m2, eff_m2(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','s','MarkerSize',8,'MarkerFaceColor',color)
		end
		axis([min(b1s_m2(1),b1s_m1(1))-.05 max(b1s_m2(end),b1s_m1(end))+.05 m-1.0 m+.2])
		grid on
		hold off
		h_leg = legend(...
							'$n_2 = 60$',...
							'$n_2 = 40$',...
							'$n_2 = 20$',...
							'$n_2 = 10$',...
							'$n_2 = 06$',...
							'Location','SouthEast');
		set(h_leg,'Interpreter','latex');
		h_xlab = xlabel('$\gamma$','Interpreter','latex');
%		h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');

		% NewCenturySchlbk, AvantGarde, Helvetica
		set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
%		set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
		set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
		box on
%		set([h_xlab, h_leg], 'FontAngle', 'Oblique')
%		title(['sigma ' num2str(s) ' - n1 ' num2str(nsims1_m1(jj))])

		print(gcf, '-depsc2', sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
		print(gcf, '-dpng',   sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	end

%	% leave only integer labels, so that final PsNR and basic PSNR figures have
%	% same width
%	ylls = get(gca,'yticklabel');
%	for ll = 1:size(ylls,1),
%		if mod(str2num(ylls(ll,:)),1) ~= 0,
%			ylls(ll,:) = ' ' * ones(1,size(ylls,2));
%		end
%	end
%	set(gca,'yticklabel', ylls(:,1:2));

	end


	if 1,
	% ------------------
	% average basic psnr
	ebb_m1(:,:,:,sidx) = mean(bb_m1,4);
	ebb_m2(:,:,:,sidx) = mean(bb_m2,4);


	m_m1 = max(max(max(ebb_m1(:,:,:,sidx))));
	m_m2 = max(max(max(ebb_m2(:,:,:,sidx))));
	m = max(m_m2, m_m1);

	for jj = 1:length(nsims1_m1),

		if nsims1_m1(jj) == 20 | nsims1_m1(jj) == 30,
			continue;
		end

		figure(fig_i), clf,
		fig_i = fig_i + 1;

		hold on
		color_step = 1/length(nsims2_m1);
		ii = 1;
		color = [0 0 0 ];
		plot(b1s_m1, ebb_m1(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','o','MarkerSize',8,'MarkerFaceColor',color)
		plot(b1s_m2, ebb_m2(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','s','MarkerSize',8,'MarkerFaceColor',color)
		axis([min(b1s_m2(1),b1s_m1(1))-.05 max(b1s_m2(end),b1s_m1(end))+.05 m-3 m+.2])
		grid on
		hold off
		h_xlab = xlabel('$\gamma$','Interpreter','latex');
%		h_ylab = ylabel('basic PSNR (dB)','Interpreter','latex');

		% NewCenturySchlbk, AvantGarde, Helvetica
		set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
%		set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
%		set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
%		set([h_xlab, h_leg], 'FontAngle', 'Oblique')
		box on

		print(gcf, '-depsc2', sprintf('%s-vs-%s_bpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
		print(gcf, '-dpng',   sprintf('%s-vs-%s_bpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	end
	end



end



