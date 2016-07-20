root_folder =  '../results/vnlbayes/';
root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/';

%stage1_folder_pos = [root_folder 'tip_table_1r_vs_1np_derf/'];
%stage1_folder_neg = [root_folder 'tip_table_1r_vs_1np_derf-neg/'];

stage1_folder_pos = [root_folder 'tip_table_1r_vs_1np_mono/'];
stage1_folder_neg = [root_folder 'tip_table_1r_vs_1np_mono-neg/'];


% stage 1 ---------------------------------------------------
% ranks_p = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,      100];
% ranks_n = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,95,100];
ranks_n = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,110,120,130,140,150,160,170,180,190,195,196];
ranks_p = [5,10,15,20,25,30,35,40,   50,   60,   70,80,   100,        130,        160,                196];
nsims = [100:100:800];
seqs = {'bus', 'football', 'foreman', 'tennis'};

eff_p = zeros(length(nsims),length(ranks_p),4);
ebb_p = zeros(length(nsims),length(ranks_p),4);

eff_n = zeros(length(nsims),length(ranks_n),4);
ebb_n = zeros(length(nsims),length(ranks_n),4);


% indices to be removed
% idx_r_n = [15,17];
% idx_r_p = [15,17,18,19];
idx_r_n = [        15,17,   19,21,   23,   25,27,   29,   31,33,   35,   37];
idx_r_p = [9,11,13,15,17,18,19,21,22,23,24,25,27,28,29,30,31,33,34,35,36,37,38,39];

fig_i = 10;
for s = 10:10:40,

	% load final neg
	ff_n = zeros(length(nsims),length(ranks_n),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_neg 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;
		tmp(1,20) = NaN; % with n = 100, r = 100 the system is ill-conditioned and we remove it

		% remove indices that have not been computed
		tmp(:,idx_r_n) = [];
		ff_n(:,:,i) = tmp;
	end

	% load basic neg
	bb_p = zeros(length(nsims),length(ranks_n),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_neg 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;
		tmp(1,20) = NaN; % with n = 100, r = 100 the system is ill-conditioned and we remove it

		tmp(:,idx_r_n) = [];
		bb_n(:,:,i) = tmp;
	end

	% load final pos
	ff_p = zeros(length(nsims),length(ranks_p),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_pos 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;
		% tmp(1,20) = NaN; % with n = 100, r = 100 the system is ill-conditioned and we remove it

		% remove indices that have not been computed
		tmp(:,idx_r_p) = [];
		ff_p(:,:,i) = tmp;

	end

	% load basic pos
	bb_p = zeros(length(nsims),length(ranks_p),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_pos 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;
		% tmp(1,20) = NaN; % with n = 100, r = 100 the system is ill-conditioned and we remove it

		tmp(:,idx_r_p) = [];
		bb_p(:,:,i) = tmp;
	end


	% ------------------
	% average final psnr
	eff_p(:,:,s/10) = mean(ff_p,3);
	eff_n(:,:,s/10) = mean(ff_n,3);

	figure(fig_i), clf,
	fig_i = fig_i + 1;

	m_p = max(max(eff_p(:,:,s/10)));
	m_n = max(max(eff_n(:,:,s/10)));
	m = max(m_n, m_p);

	hold on
	for ii = length(nsims):-1:1,
		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
		plot(ranks_p, eff_p(ii,:,s/10),'--', 'Color',color,'LineWidth', 2)
%		axis([0 105 m-3 m+.2])
		axis([0 200 m-3 m+.2])
		grid on
	end
	for ii = length(nsims):-1:1,
		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
		plot(ranks_n, eff_n(ii,:,s/10),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color)
%		axis([0 105 m-3 m+.2])
		axis([0 200 m-3 m+.2])
		grid on
	end
	hold off
	h_leg = legend('$n_1 = 800$',...
	               '$n_1 = 700$',...
	               '$n_1 = 600$',...
	               '$n_1 = 500$',...
	               '$n_1 = 400$',...
	               '$n_1 = 300$',...
	               '$n_1 = 200$',...
	               '$n_1 = 100$','Location','SouthEast');
	set(h_leg,'Interpreter','latex');
	h_xlab = xlabel('$r_1$','Interpreter','latex');
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

	print(gcf, '-depsc2', ['fpsnr_r1-np1-curves_s' num2str(s) '_average_mono_neg-vs-pos-weights']);

	% ------------------
	% average basic psnr
	ebb_p(:,:,s/10) = mean(bb_p,3);
	ebb_n(:,:,s/10) = mean(bb_n,3);

	figure(fig_i + 10), clf,
	m_p = max(max(ebb_p(:,:,s/10)));
	m_n = max(max(ebb_n(:,:,s/10)));
	m = max(m_n, m_p);

	hold on
	for ii = length(nsims):-1:1,
		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
		plot(ranks_p, ebb_p(ii,:,s/10),'--', 'Color',color,'LineWidth', 2)
%		axis([0 105 m-5 m+.2])
		axis([0 200 m-5 m+.2])
		grid on
	end
	for ii = length(nsims):-1:1,
		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
		h_plots(ii) = plot(ranks_n, ebb_n(ii,:,s/10),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color);
%		axis([0 105 m-5 m+.2])
		axis([0 200 m-5 m+.2])
		grid on
	end
	hold off
	h_leg = legend(fliplr(h_plots),...
	               '$n_1 = 800$',...
	               '$n_1 = 700$',...
	               '$n_1 = 600$',...
	               '$n_1 = 500$',...
	               '$n_1 = 400$',...
	               '$n_1 = 300$',...
	               '$n_1 = 200$',...
	               '$n_1 = 100$','Location','SouthEast');
	set(h_leg,'Interpreter','latex');
	h_xlab = xlabel('$r_1$','Interpreter','latex');
	h_ylab = ylabel('basic PSNR (dB)','Interpreter','latex');

	% NewCenturySchlbk, AvantGarde, Helvetica
	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	box on
	print(gcf, '-depsc2', ['bpsnr_r1-np1-curves_s' num2str(s) '_average_mono_neg-vs-pos-weights']);


end



