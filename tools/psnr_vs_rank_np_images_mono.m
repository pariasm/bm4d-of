root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/';

stage1_folder = [root_folder 'table_1r_vs_1np_mono/'];
stage2_folder = [root_folder 'table_2r_vs_2np_mono/'];

% % stage 2 ---------------------------------------------------
% ranks = [10:10:80];
% nsims = [50:50:200];
% nseqs = 4;
% 
% f10 = zeros(length(nsims),length(ranks),nseqs);
% f20 = zeros(length(nsims),length(ranks),nseqs);
% f40 = zeros(length(nsims),length(ranks),nseqs);
% 
% f10(:,:,1) = load([stage2_folder 'table_fpsnr_tennis_mono_s10']);
% f10(:,:,2) = load([stage2_folder 'table_fpsnr_mobile_mono_s10']);
% f10(:,:,3) = load([stage2_folder 'table_fpsnr_foreman_s10']);
% f10(:,:,4) = load([stage2_folder 'table_fpsnr_tennis_s10']);
% 
% f20(:,:,1) = load([stage2_folder 'table_fpsnr_tennis_mono_s20']);
% f20(:,:,2) = load([stage2_folder 'table_fpsnr_mobile_mono_s20']);
% f20(:,:,3) = load([stage2_folder 'table_fpsnr_foreman_s20']);
% f20(:,:,4) = load([stage2_folder 'table_fpsnr_tennis_s20']);
% 
% f40(:,:,1) = load([stage2_folder 'table_fpsnr_tennis_mono_s40']);
% f40(:,:,2) = load([stage2_folder 'table_fpsnr_mobile_mono_s40']);
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
% f10(:,:,1) = load([stage2_folder 'table_time_tennis_mono_s10']);
% f10(:,:,2) = load([stage2_folder 'table_time_mobile_mono_s10']);
% f10(:,:,3) = load([stage2_folder 'table_time_foreman_s10']);
% f10(:,:,4) = load([stage2_folder 'table_time_tennis_s10']);
% 
% f20(:,:,1) = load([stage2_folder 'table_time_tennis_mono_s20']);
% f20(:,:,2) = load([stage2_folder 'table_time_mobile_mono_s20']);
% f20(:,:,3) = load([stage2_folder 'table_time_foreman_s20']);
% f20(:,:,4) = load([stage2_folder 'table_time_tennis_s20']);
% 
% f40(:,:,1) = load([stage2_folder 'table_time_tennis_mono_s40']);
% f40(:,:,2) = load([stage2_folder 'table_time_mobile_mono_s40']);
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




%<-->% % stage 1 ---------------------------------------------------
%<-->% ranks = [10:20:210];
%<-->% nsims = [50:100:750];
%<-->% nseqs = 5;
%<-->% 
%<-->% seqs = {'tennis_mono', 'mobile_mono', 'stefan_mono', 'gflower', 'gsalesman'};
%<-->% 
%<-->% for s = 10:10:40,
%<-->% 
%<-->% 	ff = zeros(length(nsims),length(ranks),nseqs);
%<-->% 
%<-->% 	for i = 1:length(seqs),
%<-->% 		ff(:,:,i) = load([stage1_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
%<-->% 	end
%<-->% 
%<-->% 	for i = 1:length(seqs),
%<-->% 		m = max(max(ff(:,1:end-1,i)));
%<-->% 		imagesc(ff(:,:,i), m + [-.5 0]),
%<-->% 		axis equal,axis tight, colormap default
%<-->% 		colorbar
%<-->% 		set(gca,'XTick',[1:length(ranks)])
%<-->% 		set(gca,'XTickLabel',ranks)
%<-->% 		set(gca,'YTickLabel',nsims)
%<-->% 		xlabel('r')
%<-->% 		ylabel('n')
%<-->% 		print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
%<-->% 	end
%<-->% 
%<-->% 	mff = mean(ff,3);
%<-->% 	eff = mff;% - max(mff(:));
%<-->% 
%<-->% 	m = max(max(eff(:,1:end-1)));
%<-->% 	imagesc(eff, m + [-.5 0]),
%<-->% 	axis equal,axis tight, colormap default
%<-->% 	colorbar
%<-->% 	set(gca,'XTick',[1:length(ranks)])
%<-->% 	set(gca,'XTickLabel',ranks)
%<-->% 	set(gca,'YTickLabel',nsims)
%<-->% 	xlabel('r')
%<-->% 	ylabel('n')
%<-->% 	print(gcf, '-depsc2', ['fpsnr_r1-np1-image_s' num2str(s) '_average_mono']);
%<-->% 
%<-->% end
%<-->% 
%<-->% for s = 10:10:40,
%<-->% 
%<-->% 	bb = zeros(length(nsims),length(ranks),nseqs);
%<-->% 
%<-->% 	for i = 1:length(seqs),
%<-->% 		bb(:,:,i) = load([stage1_folder 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
%<-->% 	end
%<-->% 
%<-->% 	for i = 1:length(seqs),
%<-->% 		m = max(max(bb(:,1:end-1,i)));
%<-->% 		imagesc(bb(:,:,i), m + [-2 0]),
%<-->% 		axis equal,axis tight, colormap default
%<-->% 		colorbar
%<-->% 		set(gca,'XTick',[1:length(ranks)])
%<-->% 		set(gca,'XTickLabel',ranks)
%<-->% 		set(gca,'YTickLabel',nsims)
%<-->% 		xlabel('r')
%<-->% 		ylabel('n')
%<-->% 		print(gcf, '-depsc2', ['bpsnr_r1-np1-image_s' num2str(s) '_' seqs{i}]);
%<-->% 	end
%<-->% 
%<-->% 	mbb = mean(bb,3);
%<-->% 	ebb = mbb;% - max(mbb(:));
%<-->% 
%<-->% 	m = max(max(ebb(:,1:end-1)));
%<-->% 	imagesc(ebb, m + [-2 0]),
%<-->% 	axis equal,axis tight, colormap default
%<-->% 	colorbar
%<-->% 	set(gca,'XTick',[1:length(ranks)])
%<-->% 	set(gca,'XTickLabel',ranks)
%<-->% 	set(gca,'YTickLabel',nsims)
%<-->% 	xlabel('r')
%<-->% 	ylabel('n')
%<-->% 	print(gcf, '-depsc2', ['bpsnr_r1-np1-image_s' num2str(s) '_average_mono']);
%<-->% 
%<-->% end


%ranks = [10:10:100];
%nsims = [50:100:650];
%nseqs = 1;
%
%t10 = zeros(length(nsims),length(ranks),nseqs);
%t20 = zeros(length(nsims),length(ranks),nseqs);
%t30 = zeros(length(nsims),length(ranks),nseqs);
%t40 = zeros(length(nsims),length(ranks),nseqs);
%
%t10(:,:,1) = load([stage1_folder 'table_time_tennis_mono_s10']);
%t10(:,:,2) = load([stage1_folder 'table_time_mobile_mono_s10']);
%t10(:,:,3) = load([stage1_folder 'table_time_stefan_mono_s10']);
%t10(:,:,4) = load([stage1_folder 'table_time_gflower_s10']);
%t10(:,:,5) = load([stage1_folder 'table_time_gsalesman_s10']);
%
%t20(:,:,1) = load([stage1_folder 'table_time_tennis_mono_s20']);
%t20(:,:,2) = load([stage1_folder 'table_time_mobile_mono_s20']);
%t20(:,:,3) = load([stage1_folder 'table_time_stefan_mono_s20']);
%t20(:,:,4) = load([stage1_folder 'table_time_gflower_s20']);
%t20(:,:,5) = load([stage1_folder 'table_time_gsalesman_s20']);
%
%%t30(:,:,1) = load([stage1_folder 'table_time_tennis_mono_s30']);
%%t30(:,:,2) = load([stage1_folder 'table_time_mobile_mono_s30']);
%%t30(:,:,3) = load([stage1_folder 'table_time_stefan_mono_s30']);
%%t30(:,:,4) = load([stage1_folder 'table_time_gflower_s30']);
%%t30(:,:,5) = load([stage1_folder 'table_time_gsalesman_s30']);
%
%t40(:,:,1) = load([stage1_folder 'table_time_tennis_mono_s40']);
%t40(:,:,2) = load([stage1_folder 'table_time_mobile_mono_s40']);
%t40(:,:,3) = load([stage1_folder 'table_time_stefan_mono_s40']);
%t40(:,:,4) = load([stage1_folder 'table_time_gflower_s40']);
%t40(:,:,5) = load([stage1_folder 'table_time_gsalesman_s40']);
%
%mt10 = mean(t10,3);
%mt20 = mean(t20,3);
%mt30 = mean(t30,3);
%mt40 = mean(t40,3);
%
%imagesc(mt10,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s10_average_derf']);
%
%imagesc(mt20,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s20_average_derf']);
%
%imagesc(mt30,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s30_average_derf']);
%
%imagesc(mt40,[0 1000]), axis equal,axis tight, colormap gray
%set(gca,'XTickLabel',ranks)
%set(gca,'YTickLabel',nsims)
%xlabel('rank')
%ylabel('n_{sim}')
%print(gcf, '-depsc2', ['time_r1-np1-image_s40_average_derf']);



% stage 2 ---------------------------------------------------
ranks = [5:5:190];
nsims = [50:50:350];
nseqs = 5;

seqs = {'tennis_mono', 'mobile_mono', 'stefan_mono', 'gflower', 'gsalesman'};

mff = zeros(length(nsims),length(ranks),4);

ranks_l = ranks(2:4:end);
[rr_l, nn_l] = meshgrid(ranks_l,nsims);
[rr_h, nn_h] = meshgrid(ranks  ,nsims);

for s = 10,%:10:40,

	ff = zeros(length(nsims),length(ranks),nseqs);

	for i = 1:length(seqs),
		tmp = load([stage2_folder 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		tmp_l = tmp(:,2:4:end);
		ff(:,:,i) = interp2(rr_l, nn_l, tmp_l, rr_h, nn_h);
%		ff(:,:,i) = tmp(:,round((ranks - 10)/5/4)*4 + 2);
		ff(:,ranks <= 55,i) = tmp(:,ranks <= 55);
	end

	for i = 1:length(seqs),
		m = max(max(ff(:,1:end-1,i)));
		imagesc(ff(:,:,i), m + [-2 0]),
%		axis equal,
		axis tight,
		colormap default
		colorbar
		set(gca,'XTick',[1:3:length(ranks)])
		set(gca,'YTick',[1:length(nsims)])
		set(gca,'XTickLabel',ranks(1:3:end))
		set(gca,'YTickLabel',nsims)
		xlabel('r')
		ylabel('n')
		%print(gcf, '-depsc2', ['fpsnr_r2-np2-image_s' num2str(s) '_' seqs{i}]);
	end

	mff(:,:,s/10) = mean(ff,3);
	eff = mff(:,:,s/10);% - max(mff(:,:,s/10));

	m = max(max(eff(:,1:end-1)));
	imagesc(eff, m + [-2 0]),
%	axis equal,
	axis tight,
	colormap default
	colorbar
	set(gca,'XTick',[1:3:length(ranks)])
	set(gca,'YTick',[1:length(nsims)])
	set(gca,'XTickLabel',ranks(1:3:end))
	set(gca,'YTickLabel',nsims)
	xlabel('r')
	ylabel('n')
	%print(gcf, '-depsc2', ['fpsnr_r2-np2-image_s' num2str(s) '_average_mono']);

end


