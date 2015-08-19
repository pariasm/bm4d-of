root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/';

stage2_folder = [root_folder 'table_2r_vs_2np_derf/'];
stage1_folder = [root_folder 'table_1r_vs_1np_derf/'];

% % stage 2 ---------------------------------------------------
% ranks = [10:10:80];
% nsims = [50:50:200];
% nseqs = 4;
% 
% f10 = zeros(length(nsims),length(ranks),nseqs);
% f20 = zeros(length(nsims),length(ranks),nseqs);
% f40 = zeros(length(nsims),length(ranks),nseqs);
% 
% f10(:,:,1) = load([stage2_folder 'table_fpsnr_bus_s10']);
% f10(:,:,2) = load([stage2_folder 'table_fpsnr_football_s10']);
% f10(:,:,3) = load([stage2_folder 'table_fpsnr_foreman_s10']);
% f10(:,:,4) = load([stage2_folder 'table_fpsnr_tennis_s10']);
% 
% f20(:,:,1) = load([stage2_folder 'table_fpsnr_bus_s20']);
% f20(:,:,2) = load([stage2_folder 'table_fpsnr_football_s20']);
% f20(:,:,3) = load([stage2_folder 'table_fpsnr_foreman_s20']);
% f20(:,:,4) = load([stage2_folder 'table_fpsnr_tennis_s20']);
% 
% f40(:,:,1) = load([stage2_folder 'table_fpsnr_bus_s40']);
% f40(:,:,2) = load([stage2_folder 'table_fpsnr_football_s40']);
% f40(:,:,3) = load([stage2_folder 'table_fpsnr_foreman_s40']);
% f40(:,:,4) = load([stage2_folder 'table_fpsnr_tennis_s40']);
% 
% mf10 = mean(f10,3);
% mf20 = mean(f20,3);
% mf40 = mean(f40,3);
% 
% ef10 = mf10 - max(mf10(:));
% ef20 = mf20 - max(mf20(:));
% ef40 = mf40 - max(mf40(:));
% 
% imagesc(ef10,[-1 0]), axis equal,axis tight, colormap gray
% set(gca,'XTickLabel',ranks)
% set(gca,'YTickLabel',nsims)
% xlabel('rank')
% ylabel('n_{sim}')
% print(gcf, '-depsc2', ['psnr_r2-np2-image_s10_average_derf']);
% 
% imagesc(ef20,[-1 0]), axis equal,axis tight, colormap gray
% set(gca,'XTickLabel',ranks)
% set(gca,'YTickLabel',nsims)
% xlabel('rank')
% ylabel('n_{sim}')
% print(gcf, '-depsc2', ['psnr_r2-np2-image_s20_average_derf']);
% 
% imagesc(ef40,[-1 0]), axis equal,axis tight, colormap gray
% set(gca,'XTickLabel',ranks)
% set(gca,'YTickLabel',nsims)
% xlabel('rank')
% ylabel('n_{sim}')
% print(gcf, '-depsc2', ['psnr_r2-np2-image_s40_average_derf']);
% 
% 
% f10 = zeros(length(nsims),length(ranks),nseqs);
% f20 = zeros(length(nsims),length(ranks),nseqs);
% f40 = zeros(length(nsims),length(ranks),nseqs);
% 
% f10(:,:,1) = load([stage2_folder 'table_time_bus_s10']);
% f10(:,:,2) = load([stage2_folder 'table_time_football_s10']);
% f10(:,:,3) = load([stage2_folder 'table_time_foreman_s10']);
% f10(:,:,4) = load([stage2_folder 'table_time_tennis_s10']);
% 
% f20(:,:,1) = load([stage2_folder 'table_time_bus_s20']);
% f20(:,:,2) = load([stage2_folder 'table_time_football_s20']);
% f20(:,:,3) = load([stage2_folder 'table_time_foreman_s20']);
% f20(:,:,4) = load([stage2_folder 'table_time_tennis_s20']);
% 
% f40(:,:,1) = load([stage2_folder 'table_time_bus_s40']);
% f40(:,:,2) = load([stage2_folder 'table_time_football_s40']);
% f40(:,:,3) = load([stage2_folder 'table_time_foreman_s40']);
% f40(:,:,4) = load([stage2_folder 'table_time_tennis_s40']);
% 
% mf10 = mean(f10,3);
% mf20 = mean(f20,3);
% mf40 = mean(f40,3);
% 
% imagesc(mf10,[0 1000]), axis equal,axis tight, colormap gray
% set(gca,'XTickLabel',ranks)
% set(gca,'YTickLabel',nsims)
% xlabel('rank')
% ylabel('n_{sim}')
% print(gcf, '-depsc2', ['time_r2-np2-image_s10_average_derf']);
% 
% imagesc(mf20,[0 1000]), axis equal,axis tight, colormap gray
% set(gca,'XTickLabel',ranks)
% set(gca,'YTickLabel',nsims)
% xlabel('rank')
% ylabel('n_{sim}')
% print(gcf, '-depsc2', ['time_r2-np2-image_s20_average_derf']);
% 
% imagesc(mf40,[0 1000]), axis equal,axis tight, colormap gray
% set(gca,'XTickLabel',ranks)
% set(gca,'YTickLabel',nsims)
% xlabel('rank')
% ylabel('n_{sim}')
% print(gcf, '-depsc2', ['time_r2-np2-image_s40_average_derf']);




% stage 1 ---------------------------------------------------
ranks = [10:10:100];
nsims = [50:100:1250];
seqs = {'bus', 'football', 'foreman', 'tennis'};

for s = 10:10:40,

	ff = zeros(length(nsims),length(ranks),nseqs);

	for i = 1:length(seqs),
		ff(:,:,i) = load([stage1_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
	end

	for i = 1:length(seqs),
		m = max(max(ff(:,1:end-1,i)));
		imagesc(ff(:,:,i), m + [-.5 0]),
		axis equal,axis tight, colormap default
		colorbar
		set(gca,'XTick',[1:length(ranks)])
		set(gca,'YTick',[1:length(nsims)])
		set(gca,'XTickLabel',ranks)
		set(gca,'YTickLabel',nsims)
		xlabel('r')
		ylabel('n')
%		print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
	end

	mff = mean(ff,3);
	eff = mff;% - max(mff(:));

	m = max(max(eff(:,1:end-1)));
	imagesc(eff, m + [-.5 0]),
	axis equal,axis tight, colormap default
	colorbar
	set(gca,'XTick',[1:length(ranks)])
	set(gca,'YTick',[1:length(nsims)])
	set(gca,'XTickLabel',ranks)
	set(gca,'YTickLabel',nsims)
	xlabel('r')
	ylabel('n')
%	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_average_derf']);

	mff = ff - repmat(max(max(ff(:,1:end-1,:))),[length(nsims) length(ranks) 1]);
	eff = min(mff,[],3);

	m = max(max(eff(:,1:end-1)));
	imagesc(eff, m + [-.5 0]),
	axis equal,axis tight, colormap default
	colorbar
	set(gca,'XTick',[1:length(ranks)])
	set(gca,'YTick',[1:length(nsims)])
	set(gca,'XTickLabel',ranks)
	set(gca,'YTickLabel',nsims)
	xlabel('r')
	ylabel('n')
	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_min_derf']);

end

for s = 10:10:40,

	bb = zeros(length(nsims),length(ranks),nseqs);

	for i = 1:length(seqs),
		bb(:,:,i) = load([stage1_folder 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
	end

	for i = 1:length(seqs),
		m = max(max(bb(:,1:end-1,i)));
		imagesc(bb(:,:,i), m + [-2 0]),
		axis equal,axis tight, colormap default
		colorbar
		set(gca,'XTick',[1:length(ranks)])
		set(gca,'YTick',[1:length(nsims)])
		set(gca,'XTickLabel',ranks)
		set(gca,'YTickLabel',nsims)
		xlabel('r')
		ylabel('n')
%		print(gcf, '-depsc2', ['bpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
	end

	mbb = mean(bb,3);
	ebb = mbb;% - max(mbb(:));

	m = max(max(ebb(:,1:end-1)));
	imagesc(ebb, m + [-2 0]),
	axis equal,axis tight, colormap default
	colorbar
	set(gca,'XTick',[1:length(ranks)])
	set(gca,'YTick',[1:length(nsims)])
	set(gca,'XTickLabel',ranks)
	set(gca,'YTickLabel',nsims)
	xlabel('r')
	ylabel('n')
%	print(gcf, '-depsc2', ['bpsnr_r1-np1-image_s' num2str(s) '_average_derf']);

end


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


