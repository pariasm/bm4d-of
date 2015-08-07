% load data

rank = '64';
basedir2 = ['/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/'...
            'table_sizes2_1ps5x5x4_2pt4_2r' rank '_2np160/'];
basedir1 = ['/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d/results/vnlbayes/'...
            'table_sizes1_1ps5x5x4_2pt4_2r' rank '_2np160/'];
searchx = [ 13, 25, 37, 49];
patchx  = [  3,  5,  7,  9];

ylims = [33,41;  % sigma 10
         30,38;  % sigma 20
         28,36]; % sigma 40

ylimstime = [0,2000;  % sigma 10
             0,2000;  % sigma 20
             0,2000]; % sigma 40

% PSNR point cloud plots
seqs = {'army', 'dogdance', 'evergreen', 'mequon', 'walking'};
sigmas = {'10', '20', '40'};

step = 2;
plot_average = true;
plot_sequences = true;

mean_final = zeros(4,4,length(sigmas));
mean_basic = zeros(4,4,length(sigmas));
mean_times = zeros(4,4,length(sigmas));
if step == 1,
	mean_final_ref = zeros(1,length(sigmas));
	mean_basic_ref = zeros(1,length(sigmas));
	mean_times_ref = zeros(1,length(sigmas));
end

for i = 1:length(seqs),
for p = 1:length(sigmas),

	if step == 1, basedir = basedir1;
	else          basedir = basedir2;
	end

	% plots of psnr
	final = load([basedir 'table_fpsnr_' seqs{i} '_s' sigmas{p}])';
	basic = load([basedir 'table_bpsnr_' seqs{i} '_s' sigmas{p}])';
	times = load([basedir 'table_time_' seqs{i} '_s' sigmas{p}])';

	mean_final(:,:,p) = mean_final(:,:,p) + final/length(seqs);
	mean_basic(:,:,p) = mean_basic(:,:,p) + basic/length(seqs);
	mean_times(:,:,p) = mean_times(:,:,p) + times/length(seqs);


	% load reference values using default values for step 1
	if step == 1,
		basic2 = load([basedir2 'table_bpsnr_' seqs{i} '_s' sigmas{p}])';
		final2 = load([basedir2 'table_fpsnr_' seqs{i} '_s' sigmas{p}])';
		rtime2 = load([basedir2 'table_time_' seqs{i} '_s' sigmas{p}])';
		basic_ref = basic2(4,4);
		final_ref = final2(4,4);
		rtime_ref = rtime2(4,4);
		clear basic2 final2 rtime2

		mean_final_ref(p) = mean_final_ref(p) + final_ref/length(seqs);
		mean_basic_ref(p) = mean_basic_ref(p) + basic_ref/length(seqs);
		mean_times_ref(p) = mean_times_ref(p) + rtime_ref/length(seqs);
	end

	% plots of psnrs
	if plot_sequences,
		f = figure(p);
		set(f,'WindowStyle','docked');
%		bar(final')
%		set(gca, 'XTick', [1;2;3;4])
%		set(gca, 'XTickLabel', patchx)
%		legend(['r = ' num2str(searchx(1))],...
%		       ['r = ' num2str(searchx(2))],...
%		       ['r = ' num2str(searchx(3))],...
%		       ['r = ' num2str(searchx(4))],...
%				 'Location', 'Southwest');
%		xlabel('s_x')
		bar(final,1), hold on
		bar(basic, .6), 
		if step == 1,
			plot(get(gca,'XLim'),final_ref*[1,1],'m--')
			plot(get(gca,'XLim'),basic_ref*[1,1],'c--')
		end
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
		set(gca, 'XTickLabel', searchx)
		[legend_h, object_h, plot_h, text_h] = ...
		legend(['s_x = ' num2str(patchx(1))],...
				 ['s_x = ' num2str(patchx(2))],...
				 ['s_x = ' num2str(patchx(3))],...
				 ['s_x = ' num2str(patchx(4))],...
				 'Location', 'North',...
				 'Orientation', 'horizontal');
		xlabel('w_x')
		ylabel('PSNR')
		ylim(ylims(p,:))
		title(['sigma ' sigmas{p} ' ' seqs{i}])
		grid on
		box on


		if step == 1, print(gcf, '-depsc2', ['psnr_wx1-px1-bars_r' rank '_s' sigmas{p} '_' seqs{i} ])
		else          print(gcf, '-depsc2', ['psnr_wx2-px2-bars_r' rank '_s' sigmas{p} '_' seqs{i} ])
		end



		% plots of times
		f = figure(p + length(sigmas));
		set(f,'WindowStyle','docked');
%		bar(times')
%		set(gca, 'XTick', [1;2;3;4])
%		set(gca, 'XTickLabel', patchx)
%		legend(['r = ' num2str(searchx(1))],...
%		       ['r = ' num2str(searchx(2))],...
%		       ['r = ' num2str(searchx(3))],...
%		       ['r = ' num2str(searchx(4))],...
%				 'Location', 'Northwest');
%		xlabel('s_x')
		bar_h = bar(times,1);
		if step == 1,
			hold on
			plot(get(gca,'XLim'), rtime_ref*[1 1],'--')
			hold off
		end

		bars_h = findobj(gca,'Type','hggroup');
		for ii = 1:size(basic,2),
			final_color = [1-0.14*ii 0        1-0.14*ii];
			set(bars_h(ii), 'EdgeColor', final_color)
			set(bars_h(ii), 'FaceColor', final_color)
		end

		set(gca, 'XTick', [1;2;3;4])
		set(gca, 'XTickLabel', searchx)
		legend(bar_h, 4, ['s_x = ' num2str(patchx(1))],...
		                 ['s_x = ' num2str(patchx(2))],...
		                 ['s_x = ' num2str(patchx(3))],...
		                 ['s_x = ' num2str(patchx(4))],...
					  'Location', 'north',...
					  'Orientation', 'horizontal');
%		legend(['s_x = ' num2str(patchx(1))],...
%		       ['s_x = ' num2str(patchx(2))],...
%		       ['s_x = ' num2str(patchx(3))],...
%		       ['s_x = ' num2str(patchx(4))],...
%				 'Location', 'Northwest');
		xlabel('w_x')
		ylim(ylimstime(p,:))
		ylabel('computation time')
		title(['sigma ' sigmas{p} ' ' seqs{i}])
		grid on
		box on

		if step == 1, print(gcf, '-depsc2', ['time_wx1-px1-bars_r' rank '_s' sigmas{p} '_' seqs{i} ])
		else          print(gcf, '-depsc2', ['time_wx2-px2-bars_r' rank '_s' sigmas{p} '_' seqs{i} ])
		end
	end

end
%pause
end

if plot_average,

if step == 1,

	ylims = [34,39;  % sigma 10
				31,36;  % sigma 20
				28,33]; % sigma 40

	ylimstime = [0,2000;  % sigma 10
					 0,2000;  % sigma 20
					 0,2000]; % sigma 40
else

	ylims = [34,39;  % sigma 10
				32,37;  % sigma 20
				29,34]; % sigma 40

	ylimstime = [0,2000;  % sigma 10
					 0,2000;  % sigma 20
					 0,2000]; % sigma 40
end

for p = 1:length(sigmas),

	final = mean_final(:,:,p);
	basic = mean_basic(:,:,p);
	times = mean_times(:,:,p);

	% plots of psnrs
	f = figure(p);
	bar(final,1), hold on
	bar(basic,.6), 
	if step == 1,
		plot(get(gca,'XLim'),mean_final_ref(p)*[1,1],'m--')
		plot(get(gca,'XLim'),mean_basic_ref(p)*[1,1],'c--')
	end
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
	set(gca, 'XTickLabel', searchx)

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
	legend(['s_x = ' num2str(patchx(1))],...
	       ['s_x = ' num2str(patchx(2))],...
	       ['s_x = ' num2str(patchx(3))],...
	       ['s_x = ' num2str(patchx(4))],...
			 'Location', 'North',...
	       'Orientation', 'horizontal');
	xlabel('w_x')
	ylabel('PSNR')

	if step == 1, print(gcf, '-depsc2', ['psnr_wx1-px1-bars_r' rank '_s' sigmas{p} '_average' ])
	else          print(gcf, '-depsc2', ['psnr_wx2-px2-bars_r' rank '_s' sigmas{p} '_average' ])
	end

	% plots of times
	f = figure(p + length(sigmas));
	bar_h = bar(times,1);
	if step == 1,
		hold on
		plot(get(gca,'XLim'), mean_times_ref(p)*[1 1],'--')
		hold off
	end

	bars_h = findobj(gca,'Type','hggroup');
	for ii = 1:size(basic,2),
		final_color = [1-0.14*ii 0        1-0.14*ii];
		set(bars_h(ii), 'EdgeColor', final_color)
		set(bars_h(ii), 'FaceColor', final_color)
	end

	set(gca, 'XTick', [1;2;3;4])
	set(gca, 'XTickLabel', searchx)
	legend(bar_h, ['s_x = ' num2str(patchx(1))],...
	              ['s_x = ' num2str(patchx(2))],...
	              ['s_x = ' num2str(patchx(3))],...
	              ['s_x = ' num2str(patchx(4))],...
	              'Location', 'north',...
	              'Orientation', 'horizontal');
	xlabel('w_x')
	ylim(ylimstime(p,:))
	ylabel('computation time')
	set(gca, 'ygrid', 'on')
	box on

	if step == 1, print(gcf, '-depsc2', ['time_wx1-px1-bars_r' rank '_s' sigmas{p} '_average' ])
	else          print(gcf, '-depsc2', ['time_wx2-px2-bars_r' rank '_s' sigmas{p} '_average' ])
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


