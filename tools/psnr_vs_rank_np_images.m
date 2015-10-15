root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/';

%stage2_folder = [root_folder 'table_2r_vs_2np_derf/'];
stage2_folder = [root_folder 'tip_table_2r_vs_2np_derf/'];
stage1_folder = [root_folder 'tip_table_1r_vs_1np_derf/'];

% stage 2 ---------------------------------------------------

%-:-% % USE THIS CODE TO INTERPOLATE MISSING VALUES IN THE TABLES
%-:-% % THIS IS USEFULL FOR DISPLAYING RESULTS WITH IMAGESC
%-:-% ranks = [5:5:110];
%-:-% nsims = [50:50:350];
%-:-% seqs = {'bus', 'football', 'foreman', 'tennis'};
%-:-% 
%-:-% mff = zeros(length(nsims),length(ranks),length(seqs));
%-:-% 
%-:-% idx_l = [10,12,14,18,22];
%-:-% ranks_l = ranks(idx_l);
%-:-% [rr_l, nn_l] = meshgrid(ranks_l,nsims);
%-:-% [rr_h, nn_h] = meshgrid(ranks  ,nsims);
%-:-% 
%-:-% fig_i = 1;
%-:-% for s = 10:10:40,
%-:-% 
%-:-% 	ff = zeros(length(nsims),length(ranks),length(seqs));
%-:-% 
%-:-% 	for i = 1:length(seqs),
%-:-% 		tmp = load([stage2_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
%-:-% 		tmp(tmp == 0) = NaN;
%-:-% 
%-:-% 		tmp_l = tmp(:,idx_l);
%-:-% 		ff(:,:,i) = interp2(rr_l, nn_l, tmp_l, rr_h, nn_h);
%-:-% 		ff(:,ranks <= 50,i) = tmp(:,ranks <= 50);
%-:-% 	end
%-:-% 
%-:-% 	% ----------------------------
%-:-% 	% final psnr for each sequence
%-:-% 	for i = 1:length(seqs),
%-:-% 		m = max(max(ff(:,1:end-1,i)));
%-:-% 		imagesc(ff(:,:,i), m + [-.5 0]),
%-:-% 		axis equal,axis tight, colormap default
%-:-% 		colorbar
%-:-% 		set(gca,'XTick',[1:length(ranks)])
%-:-% 		set(gca,'YTick',[1:length(nsims)])
%-:-% 		set(gca,'XTickLabel',ranks)
%-:-% 		set(gca,'YTickLabel',nsims)
%-:-% 		xlabel('r')
%-:-% 		ylabel('n')
%-:-% %		print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
%-:-% 	end
%-:-% 
%-:-% 	% ------------------
%-:-% 	% average final psnr
%-:-% 	mff = mean(ff,3);
%-:-% 	eff = mff;% - max(mff(:));
%-:-% 
%-:-% 	figure(fig_i), clf, fig_i = fig_i + 1;
%-:-% 	m = max(max(eff));
%-:-% 	imagesc(eff, m + [-2 0]),
%-:-% 	axis tight, colormap default
%-:-% %	axis equal
%-:-% 	colorbar
%-:-% 	set(gca,'XTick',[1:length(ranks)])
%-:-% 	set(gca,'YTick',[1:length(nsims)])
%-:-% 	set(gca,'XTickLabel',ranks)
%-:-% 	set(gca,'YTickLabel',nsims)
%-:-% 	xlabel('r')
%-:-% 	ylabel('n')
%-:-% %	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_average_derf']);
%-:-% 
%-:-% 
%-:-% 	% ----------------------------------------
%-:-% 	% final psnr for worst performing sequence
%-:-% 	mff = ff - repmat(max(max(ff(:,1:end-1,:))),[length(nsims) length(ranks) 1]);
%-:-% 	minff = min(mff,[],3);
%-:-% 
%-:-% 	figure(fig_i), fig_i = fig_i + 1;
%-:-% 	m = max(max(minff(:,1:end-1)));
%-:-% 	imagesc(minff, m + [-2 0]),
%-:-% 	axis tight, colormap default
%-:-% %	axis equal
%-:-% 	colorbar
%-:-% 	set(gca,'XTick',[1:length(ranks)])
%-:-% 	set(gca,'YTick',[1:length(nsims)])
%-:-% 	set(gca,'XTickLabel',ranks)
%-:-% 	set(gca,'YTickLabel',nsims)
%-:-% 	xlabel('r')
%-:-% 	ylabel('n')
%-:-% %	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_min_derf']);
%-:-% 
%-:-% end

%<+>% % USE THIS CODE TO REMOVE MISSING VALUES IN THE TABLES
%<+>% % THIS IS USEFULL FOR DISPLAYING RESULTS WITH PLOT
%<+>% ranks = [5,10,15,20,25,30,35,40,45,50,60,70,90,110];
%<+>% nsims = [50:50:350];
%<+>% seqs = {'bus', 'football', 'foreman', 'tennis'};
%<+>% 
%<+>% mff = zeros(length(nsims),length(ranks),length(seqs));
%<+>% idx_r = [11,13,15,16,17,19,20,21];
%<+>% 
%<+>% fig_i = 1;
%<+>% for s = 10:10:40,
%<+>% 
%<+>% 	ff = zeros(length(nsims),length(ranks),length(seqs));
%<+>% 
%<+>% 	for i = 1:length(seqs),
%<+>% 		tmp = load([stage2_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
%<+>% 		tmp(tmp == 0) = NaN;
%<+>% 
%<+>% 		% remove indices that have not been computed
%<+>% 		tmp(:,idx_r) = [];
%<+>% 		ff(:,:,i) = tmp;
%<+>% 	end
%<+>% 
%<+>% 	% ------------------
%<+>% 	% average final psnr
%<+>% 	mff = mean(ff,3);
%<+>% 	eff = mff;% - max(mff(:));
%<+>% 
%<+>% 	figure(fig_i), clf, fig_i = fig_i + 1;
%<+>% 	m = max(max(eff));
%<+>% 
%<+>% 	hold on
%<+>% 	for ii = length(nsims):-1:1,
%<+>% 		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
%<+>% 		plot(ranks, eff(ii,:),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color)
%<+>% 		axis([0 110 m-1 m+.2])
%<+>% 		grid on
%<+>% 	end
%<+>% 	hold off
%<+>% 	h_leg = legend('n = 350',...
%<+>% 	               'n = 300',...
%<+>% 	               'n = 250',...
%<+>% 	               'n = 200',...
%<+>% 	               'n = 150',...
%<+>% 	               'n = 100',...
%<+>% 	               'n =  50','Location','SouthEast');
%<+>% 	h_xlab = xlabel('r');
%<+>% 	h_ylab = ylabel('final PSNR (dB)');
%<+>% 
%<+>% 	% NewCenturySchlbk, AvantGarde, Helvetica
%<+>% 	set([gca, h_xlab, h_ylab, h_leg], 'FontName', 'AvantGarde', 'FontSize',20)
%<+>% 	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
%<+>% 	box on
%<+>% 	print(gcf, '-depsc2', ['fpsnr_r2-np2-curves_s' num2str(s) '_average_derf']);
%<+>% 
%<+>% end

%<>% % plot time curves ---------------------------------------------------
%<>% fig_i = 1;
%<>% for s = 10:10:40,
%<>% 
%<>% 	ff = zeros(length(nsims),length(ranks),length(seqs));
%<>% 
%<>% 	for i = 1:length(seqs),
%<>% 		tmp = load([stage2_folder 'table_time_' seqs{i} '_s' num2str(s)]);
%<>% 		tmp(tmp == 0) = NaN;
%<>% 
%<>% 		tmp_l = tmp(:,idx_l);
%<>% 		ff(:,:,i) = interp2(rr_l, nn_l, tmp_l, rr_h, nn_h);
%<>% 		ff(:,ranks <= 50,i) = tmp(:,ranks <= 50);
%<>% 	end
%<>% 
%<>% 
%<>% 	% ------------------
%<>% 	% average final psnr
%<>% 	mff = median(ff,3);
%<>% 	eff = mff;% - max(mff(:));
%<>% 
%<>% 	figure(fig_i), clf, fig_i = fig_i + 1;
%<>% 	m = max(max(eff));
%<>% 
%<>% 	hold on
%<>% 	for ii = length(nsims):-1:1,
%<>% 		color = [1-0.14*(ii-1) 0 1-0.14*(ii-1)];
%<>% 		plot(ranks, eff(ii,:),'-', 'Color',color,'LineWidth', 2,'Marker','o','MarkerSize',4,'MarkerFaceColor',color)
%<>% 		axis([0 110 0 800])
%<>% 		grid on
%<>% 	end
%<>% 	hold off
%<>% 	h_leg = legend('n = 350',...
%<>% 	               'n = 300',...
%<>% 	               'n = 250',...
%<>% 	               'n = 200',...
%<>% 	               'n = 150',...
%<>% 	               'n = 100',...
%<>% 	               'n =  50','Location','SouthEast');
%<>% 	h_xlab = xlabel('r');
%<>% 	h_ylab = ylabel('final PSNR (dB)');
%<>% 
%<>% 	% NewCenturySchlbk, AvantGarde, Helvetica
%<>% 	set([gca, h_xlab, h_ylab, h_leg], 'FontName', 'AvantGarde', 'FontSize',20)
%<>% 	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
%<>% 	box on
%<>% 	print(gcf, '-depsc2', ['times_r2-np2-curves_s' num2str(s) '_average_derf']);
%<>% 
%<>% end


% stage 1 ---------------------------------------------------
ranks = [5,10,15,20,25,30,35,40,45,50,60,70,80,90,100];
nsims = [100:100:800];
seqs = {'bus', 'football', 'foreman', 'tennis'};

eff = zeros(length(nsims),length(ranks),4);
ebb = zeros(length(nsims),length(ranks),4);

% indices to be removed
idx_r = [11,13,15,17,19];

fig_i = 10;
for s = 10:10:40,

	ff = zeros(length(nsims),length(ranks),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;
		tmp(1,20) = NaN; % with n = 100, r = 100 the system is ill-conditioned and we remove it

		% remove last two ranks, 90 and 100, only for visualization
		tmp(:,19) = NaN;
		tmp(:,20) = NaN;

		% remove indices that have not been computed
		tmp(:,idx_r) = [];
		ff(:,:,i) = tmp;
	end

	bb = zeros(length(nsims),length(ranks),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;
		tmp(1,20) = NaN; % with n = 100, r = 100 the system is ill-conditioned and we remove it

		% remove last two ranks, 90 and 100, only for visualization
		tmp(:,19) = NaN;
		tmp(:,20) = NaN;

		tmp(:,idx_r) = [];
		bb(:,:,i) = tmp;
	end


	% ------------------
	% average final psnr
	mff = mean(ff,3);
	eff(:,:,s/10) = mff;% - max(mff(:));

	mbb = mean(bb,3);
	ebb(:,:,s/10) = mbb;% - max(mbb(:));

	figure(fig_i), clf, fig_i = fig_i + 1;
	m = max(max(eff(:,:,s/10)));

	hold on
	for ii = length(nsims):-1:1,
		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
		plot(ranks, eff(ii,:,s/10),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color)
		axis([0 80 m-3 m+.2])
		grid on
	end
	hold off
	h_leg = legend('n = 800',...
	               'n = 700',...
	               'n = 600',...
	               'n = 500',...
	               'n = 400',...
	               'n = 300',...
	               'n = 200',...
	               'n = 100','Location','SouthEast');
	h_xlab = xlabel('r');
	h_ylab = ylabel('final PSNR (dB)');

	% NewCenturySchlbk, AvantGarde, Helvetica
	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	box on
	print(gcf, '-depsc2', ['fpsnr_r1-np1-curves_s' num2str(s) '_average_derf']);

	figure(fig_i + 10), clf,
	m = max(max(ebb(:,:,s/10)));

	hold on
	for ii = length(nsims):-1:1,
		color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];
		plot(ranks, ebb(ii,:,s/10),'-', 'Color',color,'LineWidth', 4,'Marker','o','MarkerSize',4,'MarkerFaceColor',color)
		axis([0 80 m-5 m+.2])
		grid on
	end
	hold off
	h_leg = legend('n = 800',...
	               'n = 700',...
	               'n = 600',...
	               'n = 500',...
	               'n = 400',...
	               'n = 300',...
	               'n = 200',...
	               'n = 100','Location','SouthEast');
	h_xlab = xlabel('r');
	h_ylab = ylabel('basic PSNR (dB)');

	% NewCenturySchlbk, AvantGarde, Helvetica
	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	box on
	print(gcf, '-depsc2', ['bpsnr_r1-np1-curves_s' num2str(s) '_average_derf']);


	figure(fig_i + 20), clf,
	mb = max(max(ebb(:,:,s/10)));
	mf = max(max(eff(:,:,s/10)));

	hold on
	for jj = 1:length(ranks)-1,
		lw = (1*(jj-1) + 4*(length(ranks) - jj))/(length(ranks) - 1);
		for ii = length(nsims):-1:1,
			color = [0.14*(ii-1) 1-0.14*(ii-1) 1-0.14*(ii-1) ];

			plot(ebb(ii,jj:jj+1,s/10),eff(ii,jj:jj+1,s/10),'-', 'Color',color,'LineWidth', lw,'Marker','o','MarkerSize',lw,'MarkerFaceColor',color)
			axis([mb-.6 mb+.2 mf-.6 mf+.2])

		end

		if jj == 1,
			h_leg = legend('n = 800',...
								'n = 700',...
								'n = 600',...
								'n = 500',...
								'n = 400',...
								'n = 300',...
								'n = 200',...
								'n = 100','Location','SouthEast');
		end
	end
	grid on
	hold off
	h_xlab = xlabel('r');
	h_ylab = ylabel('basic PSNR (dB)');

	% NewCenturySchlbk, AvantGarde, Helvetica
	set([gca, h_xlab, h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	box on
	print(gcf, '-depsc2', ['bpsnr-fpsnr_r1-np1-curves_s' num2str(s) '_average_derf']);

end



% % stage 1 ---------------------------------------------------
% ranks = [10:10:100];
% nsims = [50:100:1250];
% seqs = {'bus', 'football', 'foreman', 'tennis'};
% 
% for s = 10:10:40,
% 
% 	ff = zeros(length(nsims),length(ranks),nseqs);
% 
% 	for i = 1:length(seqs),
% 		ff(:,:,i) = load([stage1_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
% 	end
% 
% 	for i = 1:length(seqs),
% 		m = max(max(ff(:,1:end-1,i)));
% 		imagesc(ff(:,:,i), m + [-.5 0]),
% 		axis equal,axis tight, colormap default
% 		colorbar
% 		set(gca,'XTick',[1:length(ranks)])
% 		set(gca,'YTick',[1:length(nsims)])
% 		set(gca,'XTickLabel',ranks)
% 		set(gca,'YTickLabel',nsims)
% 		xlabel('r')
% 		ylabel('n')
% %		print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
% 	end
% 
% 	mff = mean(ff,3);
% 	eff = mff;% - max(mff(:));
% 
% 	m = max(max(eff(:,1:end-1)));
% 	imagesc(eff, m + [-.5 0]),
% 	axis equal,axis tight, colormap default
% 	colorbar
% 	set(gca,'XTick',[1:length(ranks)])
% 	set(gca,'YTick',[1:length(nsims)])
% 	set(gca,'XTickLabel',ranks)
% 	set(gca,'YTickLabel',nsims)
% 	xlabel('r')
% 	ylabel('n')
% %	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_average_derf']);
% 
% 	mff = ff - repmat(max(max(ff(:,1:end-1,:))),[length(nsims) length(ranks) 1]);
% 	eff = min(mff,[],3);
% 
% 	m = max(max(eff(:,1:end-1)));
% 	imagesc(eff, m + [-.5 0]),
% 	axis equal,axis tight, colormap default
% 	colorbar
% 	set(gca,'XTick',[1:length(ranks)])
% 	set(gca,'YTick',[1:length(nsims)])
% 	set(gca,'XTickLabel',ranks)
% 	set(gca,'YTickLabel',nsims)
% 	xlabel('r')
% 	ylabel('n')
% 	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_min_derf']);
% 
% end
% 
% for s = 10:10:40,
% 
% 	bb = zeros(length(nsims),length(ranks),nseqs);
% 
% 	for i = 1:length(seqs),
% 		bb(:,:,i) = load([stage1_folder 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
% 	end
% 
% 	for i = 1:length(seqs),
% 		m = max(max(bb(:,1:end-1,i)));
% 		imagesc(bb(:,:,i), m + [-2 0]),
% 		axis equal,axis tight, colormap default
% 		colorbar
% 		set(gca,'XTick',[1:length(ranks)])
% 		set(gca,'YTick',[1:length(nsims)])
% 		set(gca,'XTickLabel',ranks)
% 		set(gca,'YTickLabel',nsims)
% 		xlabel('r')
% 		ylabel('n')
% %		print(gcf, '-depsc2', ['bpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
% 	end
% 
% 	mbb = mean(bb,3);
% 	ebb = mbb;% - max(mbb(:));
% 
% 	m = max(max(ebb(:,1:end-1)));
% 	imagesc(ebb, m + [-2 0]),
% 	axis equal,axis tight, colormap default
% 	colorbar
% 	set(gca,'XTick',[1:length(ranks)])
% 	set(gca,'YTick',[1:length(nsims)])
% 	set(gca,'XTickLabel',ranks)
% 	set(gca,'YTickLabel',nsims)
% 	xlabel('r')
% 	ylabel('n')
% %	print(gcf, '-depsc2', ['bpsnr_r1-np1-image_s' num2str(s) '_average_derf']);
% 
% end


%ranks = [10:10:100];
%nsims = [50:100:650];
%nseqs = 1;
%
%f10 = zeros(length(nsims),length(ranks),nseqs);
%f20 = zeros(length(nsims),length(ranks),nseqs);
%f30 = zeros(length(nsims),length(ranks),nseqs);
%f40 = zeros(length(nsims),length(ranks),nseqs);
%
%f10(:,:,1) = load([stage1_folder 'table_time_bus_s10']);
%f10(:,:,2) = load([stage1_folder 'table_time_football_s10']);
%f10(:,:,3) = load([stage1_folder 'table_time_foreman_s10']);
%f10(:,:,4) = load([stage1_folder 'table_time_tennis_s10']);
%
%f20(:,:,1) = load([stage1_folder 'table_time_bus_s20']);
%f20(:,:,2) = load([stage1_folder 'table_time_football_s20']);
%f20(:,:,3) = load([stage1_folder 'table_time_foreman_s20']);
%f20(:,:,4) = load([stage1_folder 'table_time_tennis_s20']);
%
%f30(:,:,1) = load([stage1_folder 'table_time_bus_s30']);
%f30(:,:,2) = load([stage1_folder 'table_time_football_s30']);
%f30(:,:,3) = load([stage1_folder 'table_time_foreman_s30']);
%f30(:,:,4) = load([stage1_folder 'table_time_tennis_s30']);
%
%f40(:,:,1) = load([stage1_folder 'table_time_bus_s40']);
%f40(:,:,2) = load([stage1_folder 'table_time_football_s40']);
%f40(:,:,3) = load([stage1_folder 'table_time_foreman_s40']);
%f40(:,:,4) = load([stage1_folder 'table_time_tennis_s40']);
%
%mf10 = mean(f10,3);
%mf20 = mean(f20,3);
%mf30 = mean(f30,3);
%mf40 = mean(f40,3);
%
%imagesc(mf10,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s10_average_derf']);
%
%imagesc(mf20,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s20_average_derf']);
%
%imagesc(mf30,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s30_average_derf']);
%
%imagesc(mf40,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s40_average_derf']);


