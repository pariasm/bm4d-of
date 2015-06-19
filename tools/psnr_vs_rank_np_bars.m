% load data

basedir = '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/table_rank2_1ps5x5x4_2ps5x5x4_2wx37/';

ranks = [  4,   8,  12,  16,  20];
nsim  = [ 40,  80, 120, 160, 375];

% PSNR point cloud plots
seqs = {'army', 'dogdance', 'evergreen', 'mequon', 'walking'};
sigmas = {'10', '20', '40'};

ylims = [33,40;  % sigma 10
         30,37;  % sigma 25
         28,35]; % sigma 40

ylimstime = [70,400;  % sigma 10
             70,400;  % sigma 25
             70,400]; % sigma 40

for i = 1:length(seqs),
for p = 1:length(sigmas),

	% plots of psnr
	% rows in final indicate fixed rank, columns fixed num of patches
	final = load([basedir 'table_psnr_' seqs{i} '_s' sigmas{p}]);

	f = figure(p); set(f,'WindowStyle','docked');
%	bar(final')
%	set(gca, 'XTick', [1;2;3;4;5])
%	set(gca, 'XTickLabel', nsim)
%	legend(['r = ' num2str(ranks(1))],...
%	       ['r = ' num2str(ranks(2))],...
%	       ['r = ' num2str(ranks(3))],...
%	       ['r = ' num2str(ranks(4))],...
%	       ['r = ' num2str(ranks(5))],...
%			 'Location', 'Southwest');
%	xlabel('n_{sim}')
	bar(final)
	set(gca, 'XTick', [1;2;3;4;5])
	set(gca, 'XTickLabel', ranks)
	legend(['n_{sim} = ' num2str(nsim(1))],...
	       ['n_{sim} = ' num2str(nsim(2))],...
	       ['n_{sim} = ' num2str(nsim(3))],...
	       ['n_{sim} = ' num2str(nsim(4))],...
	       ['n_{sim} = ' num2str(nsim(5))],...
			 'Location', 'Southwest');
	xlabel('rank')
	ylabel('final PSNR')
	ylim(ylims(p,:))
	title(['sigma ' sigmas{p} ' ' seqs{i}])
	grid on
	box on
%	print(gcf, '-depsc2', ['psnr_bars-final' '_s' sigmas{p} '_' seqs{i} ])



	% plots of psnr
	% rows in times indicate fixed rank, columns fixed num of patches
	times = load([basedir 'table_time_' seqs{i} '_s' sigmas{p}]);

	f = figure(p + length(sigmas)); set(f,'WindowStyle','docked');
%	bar(times')
%	set(gca, 'XTick', [1;2;3;4;5])
%	set(gca, 'XTickLabel', nsim)
%	legend(['r = ' num2str(ranks(1))],...
%	       ['r = ' num2str(ranks(2))],...
%	       ['r = ' num2str(ranks(3))],...
%	       ['r = ' num2str(ranks(4))],...
%	       ['r = ' num2str(ranks(5))],...
%			 'Location', 'Northwest');
%	xlabel('n_{sim}')
	bar(times)
	set(gca, 'XTick', [1;2;3;4;5])
	set(gca, 'XTickLabel', ranks)
	legend(['n_{sim} = ' num2str(nsim(1))],...
	       ['n_{sim} = ' num2str(nsim(2))],...
	       ['n_{sim} = ' num2str(nsim(3))],...
	       ['n_{sim} = ' num2str(nsim(4))],...
	       ['n_{sim} = ' num2str(nsim(5))],...
			 'Location', 'Northwest');
	xlabel('rank')
	ylim(ylimstime(p,:))
	ylabel('computation time')
	title(['sigma ' sigmas{p} ' ' seqs{i}])
	grid on
	box on
%	print(gcf, '-depsc2', ['psnr_bars-final' '_s' sigmas{p} '_' seqs{i} ])

end
pause
end


% for p = 1:length(sigmas),
% 	ave_time{p} = zeros(4,4);
% 	for i = 1:length(seqs),
% 		time = eval([seqs{i} '.time.time' sigmas{p}]);
% 		ave_time{p} = ave_time{p} + 1/length(seqs) * time;
% 	end
% 
% 	fid = fopen(['average_time_s' sigmas{p} '.tex'], 'w');
% 	for i = 1:4,
% 		fprintf(fid, '\t');
% 		for j = 1:3,
% 			fprintf(fid, '%6.2f & ', ave_time{p}(i,j));
% 		end
% 		fprintf(fid, '%6.2f \\\\\n', ave_time{p}(i,4));
% 	end
% 	fclose(fid);
% 
% end


