% load data

% basedir = '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/table_rank2_1ps5x5x4_2ps5x5x4_2wx37/';
% ranks = [  4,   8,  12,  16,  20];
% nsim  = [ 40,  80, 120, 160, 375];

basedir2 = '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/table_rank2_1ps5x5x4_2ps5x5x4_2wx37/';
basedir = '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/table_rank1_1ps5x5x4_2ps5x5x4_2wx37/';
ranks = [  4,   8,  12,  16,  20];
nsim  = [ 80, 120, 160, 300, 600];

ylims = [33,41;  % sigma 10
         30,38;  % sigma 25
         28,36]; % sigma 40

ylimstime = [70,400;  % sigma 10
             70,400;  % sigma 25
             70,400]; % sigma 40

% PSNR point cloud plots
seqs = {'army', 'dogdance', 'evergreen', 'mequon', 'walking'};
sigmas = {'10', '20', '40'};

for i = 1:length(seqs),
for p = 1:length(sigmas),

	% plots of psnr
	% rows in final indicate fixed rank, columns fixed num of patches
	final = load([basedir 'table_fpsnr_' seqs{i} '_s' sigmas{p}]);
	basic = load([basedir 'table_bpsnr_' seqs{i} '_s' sigmas{p}]);

	basic2 = load([basedir2 'table_bpsnr_' seqs{i} '_s' sigmas{p}]);
	final2 = load([basedir2 'table_fpsnr_' seqs{i} '_s' sigmas{p}]);
	rtime2 = load([basedir2 'table_time_' seqs{i} '_s' sigmas{p}]);
	basic_ref = basic2(4,4);
	final_ref = final2(4,4);
	rtime_ref = rtime2(4,4);
	clear basic2 final2 rtime2

	f = figure(p);
%	set(f,'WindowStyle','docked');
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
	bar(final,1), hold on
	bar(basic, .6), 
	plot(get(gca,'XLim'),final_ref*[1,1],'m--')
	plot(get(gca,'XLim'),basic_ref*[1,1],'c--')
	hold off

	bars_h = findobj(gca,'Type','hggroup');
	for ii = 1:size(basic,2),
		basic_color = [0        1-0.14*ii 1-0.14*ii];
		final_color = [1-0.14*ii 0        1-0.14*ii];
		set(bars_h(ii), 'EdgeColor', basic_color)
		set(bars_h(ii), 'FaceColor', basic_color)

		set(bars_h(ii+size(basic,2)), 'EdgeColor', final_color)
		set(bars_h(ii+size(basic,2)), 'FaceColor', final_color)
	end

	set(gca, 'XTick', [1;2;3;4;5])
	set(gca, 'XTickLabel', ranks)
	[legend_h, object_h, plot_h, text_h] = ...
	legend(['n_{sim} = ' num2str(nsim(1))],...
	       ['n_{sim} = ' num2str(nsim(2))],...
	       ['n_{sim} = ' num2str(nsim(3))],...
	       ['n_{sim} = ' num2str(nsim(4))],...
	       ['n_{sim} = ' num2str(nsim(5))],...
			 'Location', 'North',...
	       'Orientation', 'horizontal');
	xlabel('rank')
	ylabel('final PSNR')
	ylim(ylims(p,:))
	title(['sigma ' sigmas{p} ' ' seqs{i}])
	grid on
	box on


	print(gcf, '-depsc2', ['psnr_r1-np1-bars' '_s' sigmas{p} '_' seqs{i} ])



	% plots of psnr
	% rows in times indicate fixed rank, columns fixed num of patches
	times = load([basedir 'table_time_' seqs{i} '_s' sigmas{p}]);

	f = figure(p + length(sigmas));
%	set(f,'WindowStyle','docked');
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
	bar_h = bar(times,1);
	hold on
	plot(get(gca,'XLim'), rtime_ref*[1 1],'--')
	hold off

	bars_h = findobj(gca,'Type','hggroup');
	for ii = 1:size(basic,2),
		final_color = [1-0.14*ii 0        1-0.14*ii];
		set(bars_h(ii), 'EdgeColor', final_color)
		set(bars_h(ii), 'FaceColor', final_color)
	end

	set(gca, 'XTick', [1;2;3;4;5])
	set(gca, 'XTickLabel', ranks)
	legend(bar_h, 5, ['n_{sim} = ' num2str(nsim(1))],...
	           ['n_{sim} = ' num2str(nsim(2))],...
	           ['n_{sim} = ' num2str(nsim(3))],...
	           ['n_{sim} = ' num2str(nsim(4))],...
	           ['n_{sim} = ' num2str(nsim(5))],...
			     'Location', 'north',...
				  'Orientation', 'horizontal');
%	legend(['n_{sim} = ' num2str(nsim(1))],...
%	       ['n_{sim} = ' num2str(nsim(2))],...
%	       ['n_{sim} = ' num2str(nsim(3))],...
%	       ['n_{sim} = ' num2str(nsim(4))],...
%	       ['n_{sim} = ' num2str(nsim(5))],...
%			 'Location', 'Northwest');
	xlabel('rank')
	ylim(ylimstime(p,:))
	ylabel('computation time')
	title(['sigma ' sigmas{p} ' ' seqs{i}])
	grid on
	box on
	print(gcf, '-depsc2', ['time_r1-np1-bars' '_s' sigmas{p} '_' seqs{i} ])

end
%pause
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


