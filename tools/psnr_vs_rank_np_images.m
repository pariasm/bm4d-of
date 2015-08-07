root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/';

stage2_folder = [root_folder 'table_2r_vs_2np_derf/'];
stage1_folder = [root_folder 'table_1r_vs_1np_bus/'];

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
ranks = [10:10:120];
nsims = [50:100:650];
nseqs = 1;

f10 = zeros(length(nsims),length(ranks),nseqs);
f20 = zeros(length(nsims),length(ranks),nseqs);
f40 = zeros(length(nsims),length(ranks),nseqs);

f10(:,:,1) = load([stage1_folder 'table_fpsnr_bus_s10']);
%f10(:,:,2) = load([stage1_folder 'table_fpsnr_football_s10']);
%f10(:,:,3) = load([stage1_folder 'table_fpsnr_foreman_s10']);
%f10(:,:,4) = load([stage1_folder 'table_fpsnr_tennis_s10']);

f20(:,:,1) = load([stage1_folder 'table_fpsnr_bus_s20']);
%f20(:,:,2) = load([stage1_folder 'table_fpsnr_football_s20']);
%f20(:,:,3) = load([stage1_folder 'table_fpsnr_foreman_s20']);
%f20(:,:,4) = load([stage1_folder 'table_fpsnr_tennis_s20']);

f40(:,:,1) = load([stage1_folder 'table_fpsnr_bus_s40']);
%f40(:,:,2) = load([stage1_folder 'table_fpsnr_football_s40']);
%f40(:,:,3) = load([stage1_folder 'table_fpsnr_foreman_s40']);
%f40(:,:,4) = load([stage1_folder 'table_fpsnr_tennis_s40']);

ranks = [10:10:100];
f10 = f10(:,1:end-2,:)
f20 = f20(:,1:end-2,:)
f40 = f40(:,1:end-2,:)

mf10 = mean(f10,3);
mf20 = mean(f20,3);
mf40 = mean(f40,3);

ef10 = mf10 - max(mf10(:));
ef20 = mf20 - max(mf20(:));
ef40 = mf40 - max(mf40(:));

imagesc(ef10,[-1 0]), axis equal,axis tight, colormap gray
set(gca,'XTick',[1:length(ranks)])
set(gca,'XTickLabel',ranks)
set(gca,'YTickLabel',nsims)
xlabel('rank')
ylabel('n_{sim}')
print(gcf, '-depsc2', ['psnr_r1-np1-image_s10_bus']);

imagesc(ef20,[-1 0]), axis equal,axis tight, colormap gray
set(gca,'XTickLabel',ranks)
set(gca,'YTickLabel',nsims)
xlabel('rank')
ylabel('n_{sim}')
print(gcf, '-depsc2', ['psnr_r1-np1-image_s20_bus']);

imagesc(ef40,[-1 0]), axis equal,axis tight, colormap gray
set(gca,'XTickLabel',ranks)
set(gca,'YTickLabel',nsims)
xlabel('rank')
ylabel('n_{sim}')
print(gcf, '-depsc2', ['psnr_r1-np1-image_s40_bus']);


ranks = [10:10:120];
nsims = [50:100:650];
nseqs = 1;

f10 = zeros(length(nsims),length(ranks),nseqs);
f20 = zeros(length(nsims),length(ranks),nseqs);
f40 = zeros(length(nsims),length(ranks),nseqs);

f10(:,:,1) = load([stage1_folder 'table_time_bus_s10']);
%f10(:,:,2) = load([stage1_folder 'table_time_football_s10']);
%f10(:,:,3) = load([stage1_folder 'table_time_foreman_s10']);
%f10(:,:,4) = load([stage1_folder 'table_time_tennis_s10']);

f20(:,:,1) = load([stage1_folder 'table_time_bus_s20']);
%f20(:,:,2) = load([stage1_folder 'table_time_football_s20']);
%f20(:,:,3) = load([stage1_folder 'table_time_foreman_s20']);
%f20(:,:,4) = load([stage1_folder 'table_time_tennis_s20']);

f40(:,:,1) = load([stage1_folder 'table_time_bus_s40']);
%f40(:,:,2) = load([stage1_folder 'table_time_football_s40']);
%f40(:,:,3) = load([stage1_folder 'table_time_foreman_s40']);
%f40(:,:,4) = load([stage1_folder 'table_time_tennis_s40']);

ranks = [10:10:100];
f10 = f10(:,1:end-2,:)
f20 = f20(:,1:end-2,:)
f40 = f40(:,1:end-2,:)

mf10 = mean(f10,3);
mf20 = mean(f20,3);
mf40 = mean(f40,3);

imagesc(mf10,[0 1000]), axis equal,axis tight, colormap gray
set(gca,'XTickLabel',ranks)
set(gca,'YTickLabel',nsims)
xlabel('rank')
ylabel('n_{sim}')
print(gcf, '-depsc2', ['time_r1-np1-image_s10_bus']);

imagesc(mf20,[0 1000]), axis equal,axis tight, colormap gray
set(gca,'XTickLabel',ranks)
set(gca,'YTickLabel',nsims)
xlabel('rank')
ylabel('n_{sim}')
print(gcf, '-depsc2', ['time_r1-np1-image_s20_bus']);

imagesc(mf40,[0 1000]), axis equal,axis tight, colormap gray
set(gca,'XTickLabel',ranks)
set(gca,'YTickLabel',nsims)
xlabel('rank')
ylabel('n_{sim}')
print(gcf, '-depsc2', ['time_r1-np1-image_s40_bus']);


