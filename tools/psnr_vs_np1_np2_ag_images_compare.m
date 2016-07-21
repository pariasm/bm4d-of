root_folder =  '/media/pariasm/tera/funes/denoising/projects/video_nlbayes3d-dct/results/';
root_folder =  '/home/pariasm/Work/denoising/projects/video_nlbayes3d-dct/results/';

stage1_folder_m1 = [root_folder 'table_1np_vs_2np_vs_1ag_vs_2ag_img-color_1ctrY_1thY_1px8_2px4/'];

method1 = 'WIE1pos';

nsims1_m1 = [10 40 100 200];
nsims2_m1 = [10 20 40 60];
ag1s_m1  = [5 10 20 60 100 0];
ag2s_m1  = [5 10 20 60 100 0];
seqs = {...
'alley',...
'building1',...
'computer',...
'dice',...
'girl',...
'street1',...
'traffic',...
'trees'};

sigmas = [10,40];

eff_m1 = zeros(length(nsims1_m1),length(ag1s_m1),length(nsims2_m1),length(sigmas));
ebb_m1 = zeros(length(nsims1_m1),length(ag1s_m1),length(nsims2_m1),length(sigmas));

eff_m2 = zeros(length(nsims1_m1),length(ag1s_m1),length(nsims2_m1),length(sigmas));
ebb_m2 = zeros(length(nsims1_m1),length(ag1s_m1),length(nsims2_m1),length(sigmas));
                                                              

fig_i = 10;
for sidx = 1:length(sigmas), s = sigmas(sidx);

	% load final m1
	ff_m1 = nan(length(ag1s_m1),length(ag2s_m1),length(nsims1_m1),length(nsims2_m1),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_m1 'table_fpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		j = 0;
		for n1 = 1:length(nsims1_m1),
		for n2 = 1:length(nsims2_m1),
			ff_m1(:,:,n1,n2,i) = tmp([1:length(ag1s_m1)] + length(ag1s_m1)*j,:); j = j+1;
		end
		end
	end

	% load basic m1
	bb_m1 = nan(length(ag1s_m1),length(ag2s_m1),length(nsims1_m1),length(nsims2_m1),length(seqs));
	for i = 1:length(seqs),
		tmp = load([stage1_folder_m1 'table_bpsnr_' seqs{i} '_s' num2str(s)]);
		tmp(tmp == 0) = NaN;

		j = 0;
		for n1 = 1:length(nsims1_m1),
		for n2 = 1:length(nsims2_m1),
			bb_m1(:,:,n1,n2,i) = tmp([1:length(ag1s_m1)] + length(ag1s_m1)*j,:); j = j+1;
		end
		end
	end

	ag1s_m1(6) = 110;
	ag2s_m1(6) = 110;

	if 1,
	% ------------------
	% average final psnr
	eff_m1(:,:,:,sidx) = permute(squeeze(max(mean(ff_m1,5),[],2)),[2 1 3 4]);
	eff_m2(:,:,:,sidx) = permute(squeeze(min(mean(ff_m1,5),[],2)),[2 1 3 4]);


	M = max(max(max(eff_m1(:,:,:,sidx))));
	m = min(min(min(eff_m1(:,:,:,sidx))));

	for jj = 1:length(nsims1_m1),

		figure(fig_i), clf,
		fig_i = fig_i + 1;

		hold on
		color_step = 1/length(nsims2_m1);
		for ii = length(nsims2_m1):-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(ag1s_m1, eff_m1(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2)
		end
		for ii = length(nsims2_m1):-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(ag1s_m1, eff_m1(ii,:,jj,sidx),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
		end
		color_step = 1/length(nsims2_m1);
		for ii = length(nsims2_m1):-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(ag1s_m1, eff_m2(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','s','MarkerSize',8,'MarkerFaceColor',color)
		end
		axis([ag1s_m1(1)-.05 ag1s_m1(end)+.05 M-.4 M+.2])
		grid on
		hold off
		h_leg = legend(...
							'$n_2 = 60$',...
							'$n_2 = 40$',...
							'$n_2 = 20$',...
							'$n_2 = 10$',...
							'Location','SouthEast');
		set(h_leg,'Interpreter','latex');
		h_xlab = xlabel('$\gamma_{\textnormal{agg}}$','Interpreter','latex');
%		h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');

		% NewCenturySchlbk, AvantGarde, Helvetica
		set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
%		set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
		set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
		box on
%		set([h_xlab, h_leg], 'FontAngle', 'Oblique')
%		title(['sigma ' num2str(s) ' - n1 ' num2str(nsims1_m1(jj))])

		print(gcf, '-depsc2', sprintf('%s_fpsnr_np2-gagg-curves_1np%03d_s%02d_average',method1, nsims1_m1(jj), s));
		print(gcf, '-dpng',   sprintf('%s_fpsnr_np2-gagg-curves_1np%03d_s%02d_average',method1, nsims1_m1(jj), s));
	end

%	% leave only integer labels, so that final PsNR and basic PSNR figures have
%	% same width
%	ylls = get(gca,'yticklabel');
%	for ll = 1:size(ylls,1),
%		if mod(str2num(ylls(ll,:)),1) ~= 0,
%			ylls(ll,:) = ' ' * ones(1,size(ylls,2));
%		end
%	end
%	set(gca,'yticklabel', ylls(:,1:2));

	end


	if 0,
	% ------------------
	% average basic psnr
	ebb_m1(:,:,:,sidx) = mean(bb_m1,4);
	ebb_m2(:,:,:,sidx) = mean(bb_m2,4);


	m_m1 = max(max(max(ebb_m1(:,:,:,sidx))));
	m_m2 = max(max(max(ebb_m2(:,:,:,sidx))));
	m = max(m_m2, m_m1);

	for jj = 1:length(nsims1_m1),

		if nsims1_m1(jj) == 20 | nsims1_m1(jj) == 30,
			continue;
		end

		figure(fig_i), clf,
		fig_i = fig_i + 1;

		hold on
		color_step = 1/length(nsims2_m1);
		ii = 1;
		color = [0 0 0 ];
		plot(b1s_m1, ebb_m1(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','o','MarkerSize',8,'MarkerFaceColor',color)
		plot(b1s_m2, ebb_m2(ii,:,jj,sidx),'-', 'Color',color,'LineWidth', 2,'Marker','s','MarkerSize',8,'MarkerFaceColor',color)
		axis([min(b1s_m2(1),b1s_m1(1))-.05 max(b1s_m2(end),b1s_m1(end))+.05 m-3 m+.2])
		grid on
		hold off
		h_xlab = xlabel('$\gamma$','Interpreter','latex');
%		h_ylab = ylabel('basic PSNR (dB)','Interpreter','latex');

		% NewCenturySchlbk, AvantGarde, Helvetica
		set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
%		set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
%		set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
%		set([h_xlab, h_leg], 'FontAngle', 'Oblique')
		box on

		print(gcf, '-depsc2', sprintf('%s-vs-%s_bpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
		print(gcf, '-dpng',   sprintf('%s-vs-%s_bpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	end
	end



end



