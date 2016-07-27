results_folder_m1 =  '/home/pariasm/Work/denoising/projects/nldct/results/train_color_step2.bt_pos/';
results_folder_m2 =  '/home/pariasm/Work/denoising/projects/nldct/results/train_color_step2.bt_neg/';

b1s = [0.8, 0.9, 1.0, 1.1, 1.2];
%t1s = [0.0, 2.0, 2.2, 2.4, 2.6, 2.8];
t1s = [0.0, 2.0, 2.4, 2.8];
%b2s = [0.5 0.6 0.7 0.8, 0.9, 1.0, 1.1, 1.2];
b2s = [0.8, 0.9, 1.0, 1.1, 1.2];
t2s = [0.0, 2.0, 4.0, 6.0, 8.0, 10., 12.];

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

nb1s = length(b1s);
nb2s = length(b2s);
nt1s = length(t1s);
nt2s = length(t2s);

ff_m1 = nan(nb2s, nt2s, nb1s, nt1s, nsigmas, length(seqs));
bb_m1 = nan(nb2s, nt2s, nb1s, nt1s, nsigmas, length(seqs));
ff_m2 = nan(nb2s, nt2s, nb1s, nt1s, nsigmas, length(seqs));
bb_m2 = nan(nb2s, nt2s, nb1s, nt1s, nsigmas, length(seqs));
                                                              
% load data
for isig = 1:nsigmas, s = sigmas(isig);
for iseq = 1:length(seqs),

	f = load(sprintf('%stable_fpsnr_%s_s%02d', results_folder_m1, seqs{iseq}, s));
	b = load(sprintf('%stable_bpsnr_%s_s%02d', results_folder_m1, seqs{iseq}, s));

	j = 0;
	for ib1 = 1:nb1s,
	for it1 = 1:nt1s, j = j + 1;
		ff_m1(:,:,ib1,it1,isig,iseq) = f([1:nb2s] + nb2s*(j-1),:);
		bb_m1(:,:,ib1,it1,isig,iseq) = b([1:nb2s] + nb2s*(j-1),:);
	end
	end

%	f = load(sprintf('%stable_fpsnr_%s_s%02d', results_folder_m2, seqs{iseq}, s));
%	b = load(sprintf('%stable_bpsnr_%s_s%02d', results_folder_m2, seqs{iseq}, s));
%
%	j = 0;
%	for ib1 = 1:nb1s,
%	for it1 = 1:nt1s, j = j + 1;
%		ff_m2(:,:,ib1,it1,isig,iseq) = f([1:nb2s] + nb2s*(j-1),:);
%		bb_m2(:,:,ib1,it1,isig,iseq) = b([1:nb2s] + nb2s*(j-1),:);
%	end
%	end

end
end


% ignore positive t1s
ff_m1 = ff_m1(:,:,:,1,:,:);
bb_m1 = bb_m1(:,:,:,1,:,:);
t1s = [0.0];
nt1s = length(t1s);


% average psnr across sequences
mff_m1 = mean(ff_m1(:,:,:,:,:,:),6);
mbb_m1 = mean(bb_m1(:,:,:,:,:,:),6);

mff_m2 = mean(ff_m2(:,:,:,:,:,:),6);
mbb_m2 = mean(bb_m2(:,:,:,:,:,:),6);



% -----------------------------------------------------------------------------
% Effect of the search region size
%
% To analyze the effect of the search region, we plot the PSNR as a function of
% wx1 for different nsims2. In an outer loop, we consider different values of
% px2 and nsims1.
% -----------------------------------------------------------------------------

for isig = 1:nsigmas,
for ib1 = 1:nb1s,
for it1 = 1:nt1s,
	figure(isig*100 + it1*10 + ib1)
%	figure(ib1*10 + it1)
	clf

	hold on

%	color_step = 1/nt2s;
%	for ii = nt2s:-1:1,
%		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%		plot(b2s, mff_m1(:,ii,ib1,it1,isig),'-', 'Color',color,'LineWidth', 2)
%	end
%	for ii = nt2s:-1:1,
%		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%		plot(b2s, mff_m1(:,ii,ib1,it1,isig),'o','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
%	end

%	h_leg = legend(...
%						'$\tau_2 = 12.$',...
%						'$\tau_2 = 10.$',...
%						'$\tau_2 = 8.0$',...
%						'$\tau_2 = 6.0$',...
%						'$\tau_2 = 4.0$',...
%						'$\tau_2 = 2.0$',...
%						'$\tau_2 = 0.0$',...
%						'Location','SouthEast');
%	set(h_leg,'Interpreter','latex');
%	h_xlab = xlabel('$w_1$','Interpreter','latex');
%%	h_ylab = ylabel('final PSNR (dB)','Interpreter','latex');

%	for ii = nt2s:-1:1,
%		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%		plot(b2s, mff_m2(:,ii,ib1,it1,isig),'--', 'Color',color,'LineWidth', 2)
%	end
%	for ii = nt2s:-1:1,
%		color = [color_step*(ii-1) 1-color_step*(ii-1) 1-color_step*(ii-1) ];
%		plot(b2s, mff_m2(:,ii,ib1,it1,isig),'s','MarkerSize',8,'MarkerFaceColor',color, 'MarkerEdgeColor',color)
%	end

	maxPSNR_m1 = max(max(max(max(mff_m1(:,:,:,:,isig)))));
	maxPSNR_m2 = max(max(max(max(mff_m2(:,:,:,:,isig)))));
	maxPSNR = max(maxPSNR_m1, maxPSNR_m2);
	my = maxPSNR - 0.5; mx = b2s(1)   - .1;
	My = maxPSNR + 0.0; Mx = b2s(end) + .1;
	axis([mx Mx my My])
	grid on
	hold off

	imagesc(t2s, b2s, squeeze(mff_m1(:,:,ib1,it1,isig)),[my My])

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




