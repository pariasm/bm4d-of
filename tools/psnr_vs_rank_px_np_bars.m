% load data

basedir2 = ['/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/'...
            'cube_2r_2np_2px_2pt3_1px5x5x3_1rs_1nps/'];
basedir1 = ['/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/'...
            'cube_1r_1np_1px_1pt3_2px5x5x3_2r16_2np160/'];

step = 2;

patchx  = [  3,  5,  7,  9];
rank = [ 8 , 16 , 32 , 64 ];
nsim = [ 80, 160, 320, 640];

ylims = [26,41;  % sigma 10
         23,38;  % sigma 20
         21,36]; % sigma 40

ylimstime = [0,1500;  % sigma 10
             0,1500;  % sigma 20
             0,1500]; % sigma 40

% PSNR point cloud plots
seqs = {'army', 'dogdance', 'evergreen', 'walking'};
sigmas = {'10', '20', '40'};

plot_average = true;
plot_sequences = false;

mean_final = zeros(length(patchx),length(nsim),length(rank),length(sigmas));
mean_basic = zeros(length(patchx),length(nsim),length(rank),length(sigmas));
mean_times = zeros(length(patchx),length(nsim),length(rank),length(sigmas));
%if step == 2,
%	mean_final_ref = zeros(1,length(sigmas));
%	mean_basic_ref = zeros(1,length(sigmas));
%	mean_times_ref = zeros(1,length(sigmas));
%end

for i = 1:length(seqs),
for p = 1:length(sigmas),

	if step == 1, basedir = basedir1;
	else          basedir = basedir2;
	end

	% plots of psnr
	final = load([basedir 'cube_fpsnr_' seqs{i} '_s' sigmas{p}]);
	basic = load([basedir 'cube_bpsnr_' seqs{i} '_s' sigmas{p}]);
	times = load([basedir 'cube_time_'  seqs{i} '_s' sigmas{p}]);

	final = reshape(final', [length(nsim), length(rank), length(patchx)]);
	basic = reshape(basic', [length(nsim), length(rank), length(patchx)]);
	times = reshape(times', [length(nsim), length(rank), length(patchx)]);

	final = permute(final, [3 1 2]);
	basic = permute(basic, [3 1 2]);
	times = permute(times, [3 1 2]);
	% Dimensions at this point: px (rows) x nsim (colums) x rank (slices)

	mean_final(:,:,:,p) = mean_final(:,:,:,p) + final/length(seqs);
	mean_basic(:,:,:,p) = mean_basic(:,:,:,p) + basic/length(seqs);
	mean_times(:,:,:,p) = mean_times(:,:,:,p) + times/length(seqs);

	% load reference values using default values for step 1
	%if step == 2,
	%	basic1 = load([basedir1 'cube_bpsnr_' seqs{i} '_s' sigmas{p}])';
	%	final1 = load([basedir1 'cube_fpsnr_' seqs{i} '_s' sigmas{p}])';
	%	rtime1 = load([basedir1 'cube_time_'  seqs{i} '_s' sigmas{p}])';

	%	final1 = permute(reshape(final1', [length(nsim), length(rank), length(patchx)]), [2 1 3]);
	%	basic1 = permute(reshape(basic1', [length(nsim), length(rank), length(patchx)]), [2 1 3]);
	%	times1 = permute(reshape(times1', [length(nsim), length(rank), length(patchx)]), [2 1 3]);

	%	basic_ref = basic2(4,4);
	%	final_ref = final2(4,4);
	%	rtime_ref = rtime2(4,4);
	%	clear basic2 final2 rtime2

	%	mean_final_ref(p) = mean_final_ref(p) + final_ref/length(seqs);
	%	mean_basic_ref(p) = mean_basic_ref(p) + basic_ref/length(seqs);
	%	mean_times_ref(p) = mean_times_ref(p) + rtime_ref/length(seqs);
	%end

	% plots of psnrs
	if plot_sequences,
	for j = 1:length(rank),

		f = figure(1);
		set(f,'WindowStyle','docked');
		bar(final(:,:,j),1), hold on
		bar(basic(:,:,j), .6), 
		%if step == 2,
		%	plot(get(gca,'XLim'),final_ref*[1,1],'m--')
		%	plot(get(gca,'XLim'),basic_ref*[1,1],'c--')
		%end
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

		set(gca, 'XTick', [1;2;3;4])
		set(gca, 'XTickLabel', patchx)
		xlabel('s_x')
		ylabel('PSNR')
		ylim(ylims(p,:))
		[legend_h, object_h, plot_h, text_h] = ...
		legend(['n_{sim} = ' num2str(nsim(1))],...
				 ['n_{sim} = ' num2str(nsim(2))],...
				 ['n_{sim} = ' num2str(nsim(3))],...
				 ['n_{sim} = ' num2str(nsim(4))],...
				 'Location', 'North',...
				 'Orientation', 'horizontal');
		title(['sigma ' sigmas{p} ' ' seqs{i} ' r ' num2str(rank(j))])
		grid on
		box on

		if step == 1, print(gcf, '-depsc2', ['psnr_px1-np1-bars_1r' num2str(rank(j)) '_s' sigmas{p} '_' seqs{i} ])
		else          print(gcf, '-depsc2', ['psnr_px2-np2-bars_2r' num2str(rank(j)) '_s' sigmas{p} '_' seqs{i} ])
		end



		% plots of times
		f = figure(1 + length(sigmas));
		set(f,'WindowStyle','docked');
		bar_h = bar(times(:,:,j),1);
		%if step == 2,
		%	hold on
		%	plot(get(gca,'XLim'), rtime_ref*[1 1],'--')
		%	hold off
		%end

		bars_h = findobj(gca,'Type','hggroup');
		for ii = 1:size(basic,2),
			final_color = [1-0.14*ii 0        1-0.14*ii];
			set(bars_h(ii), 'EdgeColor', final_color)
			set(bars_h(ii), 'FaceColor', final_color)
		end

		set(gca, 'XTick', [1;2;3;4])
		set(gca, 'XTickLabel', patchx)
		xlabel('s_x')
		ylim(ylimstime(p,:))
		ylabel('computation time')
		legend(bar_h, 4, ['n_{sim} = ' num2str(nsim(1))],...
		                 ['n_{sim} = ' num2str(nsim(2))],...
		                 ['n_{sim} = ' num2str(nsim(3))],...
		                 ['n_{sim} = ' num2str(nsim(4))],...
					  'Location', 'north',...
					  'Orientation', 'horizontal');
		title(['sigma ' sigmas{p} ' ' seqs{i} ' r ' num2str(rank(j))])
		grid on
		box on

		if step == 1, print(gcf, '-depsc2', ['time_px1-np1-bars_1r' num2str(rank(j)) '_s' sigmas{p} '_' seqs{i} ])
		else          print(gcf, '-depsc2', ['time_px2-np2-bars_2r' num2str(rank(j)) '_s' sigmas{p} '_' seqs{i} ])
		end
		pause(.5)
	end
	end

end
end

if plot_average,

if step == 1,

	ylims = [33,38;  % sigma 10
				31,36;  % sigma 20
				28,33]; % sigma 40

	ylimstime = [0,1500;  % sigma 10
					 0,1500;  % sigma 20
					 0,1500]; % sigma 40
else

	ylims = [29,39;  % sigma 10
				27,37;  % sigma 20
				24,34]; % sigma 40

	ylimstime = [0,1500;  % sigma 10
					 0,1500;  % sigma 20
					 0,1500]; % sigma 40
end

for p = 1:length(sigmas),

	final = mean_final(:,:,:,p);
	basic = mean_basic(:,:,:,p);
	times = mean_times(:,:,:,p);

	% plots of psnrs
	for j = 1:length(rank),

		f = figure(1);
		bar(final(:,:,j),1), hold on
		bar(basic(:,:,j),.6), 
		%if step == 2,
		%	plot(get(gca,'XLim'),mean_final_ref(p)*[1,1],'m--')
		%	plot(get(gca,'XLim'),mean_basic_ref(p)*[1,1],'c--')
		%end
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

		xlabel('s_x')
		set(gca, 'XTick', [1;2;3;4])
		set(gca, 'XTickLabel', patchx)

		ylabel('PSNR')
		ylim(ylims(p,:))
		yticks = get(gca,'YTickLabel');
		for t = 1:size(yticks,1),
			if (mod(str2double(yticks(t,:)),1) ~= 0)
				yticks(t,:) = '    ';
			end
		end
		set(gca,'YTickLabel',yticks)

		set(gca, 'ygrid', 'on')
		box on

		[legend_h, object_h, plot_h, text_h] = ...
		legend(['n_{sim} = ' num2str(nsim(1))],...
				 ['n_{sim} = ' num2str(nsim(2))],...
				 ['n_{sim} = ' num2str(nsim(3))],...
				 ['n_{sim} = ' num2str(nsim(4))],...
				 'Location', 'North',...
				 'Orientation', 'horizontal');

		if step == 1, print(gcf, '-depsc2', ['psnr_px1-np1-bars_1r' num2str(rank(j)) '_s' sigmas{p} '_average' ])
		else          print(gcf, '-depsc2', ['psnr_px2-np2-bars_2r' num2str(rank(j)) '_s' sigmas{p} '_average' ])
		end

		% plots of times
		f = figure(1 + length(sigmas));
		bar_h = bar(times(:,:,j),1);
		%if step == 2,
		%	hold on
		%	plot(get(gca,'XLim'), mean_times_ref(p)*[1 1],'--')
		%	hold off
		%end

		bars_h = findobj(gca,'Type','hggroup');
		for ii = 1:size(basic,2),
			final_color = [1-0.14*ii 0        1-0.14*ii];
			set(bars_h(ii), 'EdgeColor', final_color)
			set(bars_h(ii), 'FaceColor', final_color)
		end

		xlabel('s_x')
		set(gca, 'XTick', [1;2;3;4])
		set(gca, 'XTickLabel', patchx)

		legend(bar_h, ['n_{sim} = ' num2str(nsim(1))],...
						  ['n_{sim} = ' num2str(nsim(2))],...
						  ['n_{sim} = ' num2str(nsim(3))],...
						  ['n_{sim} = ' num2str(nsim(4))],...
						  'Location', 'north',...
						  'Orientation', 'horizontal');
		ylabel('computation time')
		ylim(ylimstime(p,:))
		set(gca, 'ygrid', 'on')
		box on

		if step == 1, print(gcf, '-depsc2', ['time_px1-np1-bars_1r' num2str(rank(j)) '_s' sigmas{p} '_average' ])
		else          print(gcf, '-depsc2', ['time_px2-np2-bars_2r' num2str(rank(j)) '_s' sigmas{p} '_average' ])
		end
	end

end
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


