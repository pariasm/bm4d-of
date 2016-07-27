results_folder_m1 =  '/home/pariasm/Work/denoising/projects/nldct/results/train_gray_step1_pos/';
results_folder_m2 =  '/home/pariasm/Work/denoising/projects/nldct/results/train_gray_step1_neg/';

nsims1 = [20 30 40 50 60];
nsims2 = [10 20 40];
px1s = [4 6 8 12];
wx1s = [11 21 31 41];

seqs = {...
'gray_alley',...
'gray_building1',...
'gray_computer',...
'gray_dice',...
'gray_girl',...
'gray_street1',...
'gray_traffic',...
'gray_trees'};

sigmas = [5,10,20,30,40,60];

nsigmas = length(sigmas);
nnsims1 = length(nsims1);
nnsims2 = length(nsims2);
npx1s   = length(px1s);
nwx1s   = length(wx1s);

ff_m1 = nan(nwx1s,nnsims1,npx1s,nnsims2,nsigmas,length(seqs));
bb_m1 = nan(nwx1s,nnsims1,npx1s,nnsims2,nsigmas,length(seqs));
ff_m2 = nan(nwx1s,nnsims1,npx1s,nnsims2,nsigmas,length(seqs));
bb_m2 = nan(nwx1s,nnsims1,npx1s,nnsims2,nsigmas,length(seqs));
                                                              
% load data
for isig = 1:nsigmas, s = sigmas(isig);
for iseq = 1:length(seqs),
	f = load(sprintf('%stable_fpsnr_%s_s%02d', results_folder_m1, seqs{iseq}, s));
	b = load(sprintf('%stable_bpsnr_%s_s%02d', results_folder_m1, seqs{iseq}, s));

	j = 0;
	for inp2 = 1:nnsims2,
	for ipx1 = 1:npx1s, j = j + 1;
		ff_m1(:,:,ipx1,inp2,isig,iseq) = f([1:nwx1s] + nwx1s*(j-1),:);
		bb_m1(:,:,ipx1,inp2,isig,iseq) = b([1:nwx1s] + nwx1s*(j-1),:);
	end
	end

%	f = load(sprintf('%stable_fpsnr_%s_s%02d', results_folder_m2, seqs{iseq}, s));
%	b = load(sprintf('%stable_bpsnr_%s_s%02d', results_folder_m2, seqs{iseq}, s));
%
%	j = 0;
%	for inp2 = 1:nnsims2,
%	for ipx1 = 1:npx1s, j = j + 1;
%		ff_m2(:,:,ipx1,inp2,isig,iseq) = f([1:nwx1s] + nwx1s*(j-1),:);
%		bb_m2(:,:,ipx1,inp2,isig,iseq) = b([1:nwx1s] + nwx1s*(j-1),:);
%	end
%	end
end
end

if 0,
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full viz (a simpler viz is bellow)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% average psnr across sequences
	mff_m1 = mean(ff_m1,6);
	mbb_m1 = mean(bb_m1,6);
	
	mff_m2 = mean(ff_m2,6);
	mbb_m2 = mean(bb_m2,6);
	
	
	% -----------------------------------------------------------------------------
	% Effect of the search region size
	%
	% To analyze the effect of the search region, we plot the PSNR as a function of
	% wx1 for different nsims1. In an outer loop, we consider different values of
	% px1 and nsims2.
	% -----------------------------------------------------------------------------
	
	for isig = 1,%1:nsigmas,
	for inp2 = 1:nnsims2,
	for ipx1 = 1:npx1s,
	%	figure(isig*100 + inp2*10 + ipx1)
		figure(500 + inp2*10 + ipx1)
		clf
	
		hold on
	
		color_step = 1/nnsims1;
		for ii = nnsims1:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(wx1s, mff_m1(:,ii,ipx1,inp2,isig),'-', 'Color',color,'LineWidth', 2)
		end
		for ii = nnsims1:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(wx1s, mff_m1(:,ii,ipx1,inp2,isig),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
		end
	
	%	h_leg = legend(...
	%						'$n_1 = 60$',...
	%						'$n_1 = 50$',...
	%						'$n_1 = 40$',...
	%						'$n_1 = 30$',...
	%						'$n_1 = 20$',...
	%						'Location','SouthEast');
	%	set(h_leg,'Interpreter','latex');
	%	h_xlab = xlabel('$w_1$','Interpreter','latex');
	%%	h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');
	
		for ii = nnsims1:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(wx1s, mff_m2(:,ii,ipx1,inp2,isig),'--', 'Color',color,'LineWidth', 2)
		end
		for ii = nnsims1:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(wx1s, mff_m2(:,ii,ipx1,inp2,isig),'s','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
		end
	
		maxPSNR_m1 = max(max(max(max(mff_m1(:,:,:,:,isig)))));
		maxPSNR_m2 = max(max(max(max(mff_m2(:,:,:,:,isig)))));
		maxPSNR = max(maxPSNR_m1, maxPSNR_m2);
		my = maxPSNR - 1.0; mx = wx1s(1)   - 5;
		My = maxPSNR + 0.1; Mx = wx1s(end) + 5;
		axis([mx Mx my My])
		grid on
		hold off
	
		% NewCenturySchlbk, AvantGarde, Helvetica
	%	set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
	%	set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	%	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
		box on
	%	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	%	title(['sigma ' num2str(s) ' - n1 ' num2str(nsims1_m1(jj))])
	
	%	print(gcf, '-depsc2', sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	%	print(gcf, '-dpng',   sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	
	end
	end
	end
else
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Better viz : remove wx1 from the pixture
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	ff_m1 = squeeze(ff_m1(4,:,:,:,:,:));
	ff_m2 = squeeze(ff_m2(4,:,:,:,:,:));
	bb_m1 = squeeze(bb_m1(4,:,:,:,:,:));
	bb_m2 = squeeze(bb_m2(4,:,:,:,:,:));
	
	% average psnr across sequences
	mff_m1 = mean(ff_m1,5);
	mbb_m1 = mean(bb_m1,5);
	
	mff_m2 = mean(ff_m2,5);
	mbb_m2 = mean(bb_m2,5);
	
	
	% -----------------------------------------------------------------------------
	% Effect of the search region size
	%
	% To analyze the effect of the search region, we plot the PSNR as a function of
	% wx1 for different nsims1. In an outer loop, we consider different values of
	% px1 and nsims2.
	% -----------------------------------------------------------------------------
	
	for isig = 1:nsigmas,
	for inp2 = 1:nnsims2,
	%	figure(isig*100 + inp2*10 + ipx1)
		figure(100*isig + inp2)
		clf
	
		hold on
	
%		color_step = 1/nnsims1;
%		for ii = nnsims1:-1:1,
%			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%			plot(px1s, mff_m1(ii,:,inp2,isig),'-', 'Color',color,'LineWidth', 2)
%		end
%		for ii = nnsims1:-1:1,
%			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%			plot(px1s, mff_m1(ii,:,inp2,isig),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
%		end
		color_step = 1/npx1s;
		for ii = npx1s:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(nsims1, mff_m1(:,ii,inp2,isig),'-', 'Color',color,'LineWidth', 2)
		end
		for ii = npx1s:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(nsims1, mff_m1(:,ii,inp2,isig),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
		end
	
	%	h_leg = legend(...
	%						'$n_1 = 60$',...
	%						'$n_1 = 50$',...
	%						'$n_1 = 40$',...
	%						'$n_1 = 30$',...
	%						'$n_1 = 20$',...
	%						'Location','SouthEast');
	%	set(h_leg,'Interpreter','latex');
	%	h_xlab = xlabel('$w_1$','Interpreter','latex');
	%%	h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');
	
		for ii = nnsims1:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(px1s, mff_m2(ii,:,inp2,isig),'--', 'Color',color,'LineWidth', 2)
		end
		for ii = nnsims1:-1:1,
			color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
			plot(px1s, mff_m2(ii,:,inp2,isig),'s','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
		end
	
		maxPSNR_m1 = max(max(max(max(mff_m1(:,:,:,isig)))));
		maxPSNR_m2 = max(max(max(max(mff_m2(:,:,:,isig)))));
		maxPSNR = max(maxPSNR_m1, maxPSNR_m2);
%		my = maxPSNR - 1.0; mx = px1s(1)   - 2;
%		My = maxPSNR + 0.1; Mx = px1s(end) + 2;
		my = maxPSNR - 1.0; mx = nsims1(1)   - 10;
		My = maxPSNR + 0.1; Mx = nsims1(end) + 10;
		axis([mx Mx my My])
		grid on
		hold off
	
		% NewCenturySchlbk, AvantGarde, Helvetica
	%	set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
	%	set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
	%	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
		box on
	%	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
	%	title(['sigma ' num2str(s) ' - n1 ' num2str(nsims1_m1(jj))])
	
	%	print(gcf, '-depsc2', sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	%	print(gcf, '-dpng',   sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims1_m1(jj), s));
	
	end
	end
end
