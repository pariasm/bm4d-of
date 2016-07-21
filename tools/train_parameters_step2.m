results_folder_m1 =  '/home/pariasm/Work/denoising/projects/video_nlbayes3d-dct/results/train_color_step2_neg/';
results_folder_m2 =  '/home/pariasm/Work/denoising/projects/video_nlbayes3d-dct/results/train_color_step2_pos/';

nsims1 = [10, 20, 30, 40, 60];
nsims2 = [5, 10, 15, 20, 30];
px2s  = [3, 4, 6];
wx2s  = [11, 21, 31, 41];

%nsims1 = [20 30 40 50 60];
%nsims2 = [10 20 40];
%px1s = [4 6 8 12];
%wx1s = [11 21 31 41];

seqs = {...
'alley',...
'building1',...
'computer',...
'dice',...
'girl',...
'street1',...
'traffic',...
'trees'};

sigmas = [5,10,20,30,40];

nsigmas = length(sigmas);
nnsims2 = length(nsims2);
nnsims1 = length(nsims1);
npx2s   = length(px2s);
nwx2s   = length(wx2s);

ff_m1 = nan(nwx2s,nnsims2,npx2s,nnsims1,nsigmas,length(seqs));
bb_m1 = nan(nwx2s,nnsims2,npx2s,nnsims1,nsigmas,length(seqs));
ff_m2 = nan(nwx2s,nnsims2,npx2s,nnsims1,nsigmas,length(seqs));
bb_m2 = nan(nwx2s,nnsims2,npx2s,nnsims1,nsigmas,length(seqs));
                                                              
% load data
for isig = 1:nsigmas, s = sigmas(isig);
for iseq = 1:length(seqs),
	f = load(sprintf('%stable_fpsnr_%s_s%02d', results_folder_m1, seqs{iseq}, s));
	b = load(sprintf('%stable_bpsnr_%s_s%02d', results_folder_m1, seqs{iseq}, s));

	j = 0;
	for inp1 = 1:nnsims1,
	for ipx2 = 1:npx2s, j = j + 1;
		ff_m1(:,:,ipx2,inp1,isig,iseq) = f([1:nwx2s] + nwx2s*(j-1),:);
		bb_m1(:,:,ipx2,inp1,isig,iseq) = b([1:nwx2s] + nwx2s*(j-1),:);
	end
	end

%	f = load(sprintf('%stable_fpsnr_%s_s%02d', results_folder_m2, seqs{iseq}, s));
%	b = load(sprintf('%stable_bpsnr_%s_s%02d', results_folder_m2, seqs{iseq}, s));
%
%	j = 0;
%	for inp1 = 1:nnsims1,
%	for ipx2 = 1:npx2s, j = j + 1;
%		ff_m2(:,:,ipx2,inp1,isig,iseq) = f([1:nwx2s] + nwx2s*(j-1),:);
%		bb_m2(:,:,ipx2,inp1,isig,iseq) = b([1:nwx2s] + nwx2s*(j-1),:);
%	end
%	end
end
end

% average psnr across sequences
mff_m1 = mean(ff_m1,6);
mbb_m1 = mean(bb_m1,6);

mff_m2 = mean(ff_m2,6);
mbb_m2 = mean(bb_m2,6);


% -----------------------------------------------------------------------------
% Effect of the search region size
%
% To analyze the effect of the search region, we plot the PSNR as a function of
% wx1 for different nsims2. In an outer loop, we consider different values of
% px2 and nsims1.
% -----------------------------------------------------------------------------

for isig = 1,%1:nsigmas,
for inp1 = 1:nnsims1,
for ipx2 = 1:npx2s,
%	figure(isig*100 + inp1*10 + ipx2)
	figure(500 + inp1*10 + ipx2)
	clf

	hold on

	color_step = 1/nnsims2;
	for ii = nnsims2:-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(wx2s, mff_m1(:,ii,ipx2,inp1,isig),'-', 'Color',color,'LineWidth', 2)
	end
	for ii = nnsims2:-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(wx2s, mff_m1(:,ii,ipx2,inp1,isig),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
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

	for ii = nnsims2:-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(wx2s, mff_m2(:,ii,ipx2,inp1,isig),'--', 'Color',color,'LineWidth', 2)
	end
	for ii = nnsims2:-1:1,
		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
		plot(wx2s, mff_m2(:,ii,ipx2,inp1,isig),'s','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
	end

	maxPSNR_m1 = max(max(max(max(mff_m1(:,:,:,:,isig)))));
	maxPSNR_m2 = max(max(max(max(mff_m2(:,:,:,:,isig)))));
	maxPSNR = max(maxPSNR_m1, maxPSNR_m2);
	my = maxPSNR - 1.0; mx = wx2s(1)   - 5;
	My = maxPSNR + 0.1; Mx = wx2s(end) + 5;
	axis([mx Mx my My])
	grid on
	hold off

	% NewCenturySchlbk, AvantGarde, Helvetica
%	set([gca, h_xlab], 'FontName', 'AvantGarde', 'FontSize',20)
%	set([h_ylab], 'FontName', 'AvantGarde', 'FontSize',20)
%	set([h_leg], 'FontName', 'AvantGarde', 'FontSize',18)
	box on
%	set([h_xlab, h_leg], 'FontAngle', 'Oblique')
%	title(['sigma ' num2str(s) ' - n1 ' num2str(nsims2_m1(jj))])

%	print(gcf, '-depsc2', sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims2_m1(jj), s));
%	print(gcf, '-dpng',   sprintf('%s-vs-%s_fpsnr_np2-g-curves_1np%03d_s%02d_average',method1, method2, nsims2_m1(jj), s));

end
end
end




