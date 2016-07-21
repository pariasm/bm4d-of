root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d-dct/results/';
root_folder =  '/home/pariasm/Work/denoising/projects/video_nlbayes3d-dct/results/';


% stage1_folder_pos = [root_folder 'table_1np_vs_2np_img-color_1cntrY_2cntrY_1px8_2px8_pos/'];
% stage1_folder_neg = [root_folder 'table_1np_vs_2np_img-color_1cntrN_2cntrY_1px8_2px8_pos/'];

stage1_folder_pos = [root_folder 'table_1np_vs_2np_img-color_1cntrY_2cntrY_1px8_2px8_pos/'];
stage1_folder_neg = [root_folder 'table_1np_vs_2np_img-color_1cntrY_2cntrN_1px8_2px8_pos/'];


% stage 1 ---------------------------------------------------
nsims1_p = [4 6 8 10 15 20 30 40 50 60 80 100];
nsims1_n = [4 6 8 10 15 20 30 40 50 60 80 100];
nsims2_p = [4 6 8 10 15 20 30 40 50 60 80 100];
nsims2_n = [1 2 3 4 6 8 10 15 20 30 40 50 60 80 100];
seqs = {...
'alley',...
'book',...
'building1',...
'building2',...
'computer',...
'dice',...
'flowers1',...
'flowers2',...
'gardens',...
'girl',...
'hallway',...
'man1',...
'man2',...
'plaza',...
'statue',...
'street1',...
'street2',...
'traffic',...
'trees',...
'valldemossa',...
'yard'};

sigmas = [10,20,40];

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
		plot(nsims2_p, eff_p(ii,:,sidx),'--', 'Color',color,'LineWidth', 2)
		axis([0 100 m-3 m+.2])
		grid on
	end
	color_step = 1/length(nsims1_n);
	for ii = length(nsims1_n):-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(nsims2_n, eff_n(ii,:,sidx),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color)
		axis([0 100 m-3 m+.2])
		grid on
	end
	hold off
	h_leg = legend('$n_1 = 100$',...
	               '$n_1 =  80$',...
	               '$n_1 =  60$',...
	               '$n_1 =  50$',...
	               '$n_1 =  40$',...
	               '$n_1 =  30$',...
	               '$n_1 =  20$',...
	               '$n_1 =  15$',...
	               '$n_1 =  10$',...
	               '$n_1 =   8$',...
	               '$n_1 =   6$',...
	               '$n_1 =   4$','Location','SouthEast');
	set(h_leg,'Interpreter','latex');
	h_xlab = xlabel('$n_2$','Interpreter','latex');
	h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');

	% NewCenturySchlbk, AvantGarde, Helvetica
	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	box on

	% leave only integer labels, so that final PsNR and basic PSNR figures have
	% same width
	ylls = get(gca,'yticklabel');
	for ll = 1:size(ylls,1),
		if mod(str2num(ylls(ll,:)),1) ~= 0,
			ylls(ll,:) = ' ' * ones(1,size(ylls,2));
		end
	end
	set(gca,'yticklabel', ylls(:,1:2));

%	print(gcf, '-depsc2', ['fpsnr_r1-np1-curves_s' num2str(s) '_average_mono_neg-vs-pos-weights']);

	% ------------------
	% average basic psnr
	ebb_p(:,:,sidx) = mean(bb_p,3);
	ebb_n(:,:,sidx) = mean(bb_n,3);

%	figure(fig_i + 10), clf,
%	m_p = max(max(ebb_p(:,:,sidx)));
%	m_n = max(max(ebb_n(:,:,sidx)));
%	m = max(m_n, m_p);
%
%	hold on
%	color_step = 1/length(nsims1_p);
%	for ii = length(nsims1_p):-1:1,
%		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%		plot(nsims2_p, ebb_p(ii,:,sidx),'--', 'Color',color,'LineWidth', 2)
%		axis([0 100 m-5 m+.2])
%		grid on
%	end
%	color_step = 1/length(nsims1_p);
%	for ii = length(nsims1_n):-1:1,
%		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%		h_plots(ii) = plot(nsims2_n, ebb_n(ii,:,sidx),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color);
%		axis([0 100 m-5 m+.2])
%		grid on
%	end
%	hold off
%	h_leg = legend('$n_1 = 100$',...
%	               '$n_1 =  80$',...
%	               '$n_1 =  60$',...
%	               '$n_1 =  50$',...
%	               '$n_1 =  40$',...
%	               '$n_1 =  30$',...
%	               '$n_1 =  20$',...
%	               '$n_1 =  15$',...
%	               '$n_1 =  10$',...
%	               '$n_1 =   8$',...
%	               '$n_1 =   6$',...
%	               '$n_1 =   4$','Location','SouthEast');
%	set(h_leg,'Interpreter','latex');
%	h_xlab = xlabel('$n_2$','Interpreter','latex');
%	h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');
%
%	% NewCenturySchlbk, AvantGarde, Helvetica
%	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
%	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
%	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
%	box on
%%	print(gcf, '-depsc2', ['bpsnr_r1-np1-curves_s' num2str(s) '_average_mono_neg-vs-pos-weights']);


end



