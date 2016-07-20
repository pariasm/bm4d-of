% In this script we plot the pdfs of the variances and of the weights
% and their Gaussian approximation

% parameters
snr = .05;
sigma = 20;
N = 20;

%<>% % plot results varying N
%<>% maxh = 0;
%<>% Ns = 2.^[9:-2:3];
%<>% hplots = zeros(1,length(Ns));
%<>% for ii = 1:length(Ns), N = Ns(ii); 
%<>% 
%<>% 	% first, a scaled non-central chi square
%<>% 	
%<>% 	v = [0:.1:2000];
%<>% 	
%<>% 	mu = N*snr;
%<>% 	a = sigma^2/N;
%<>% 	
%<>% 	
%<>% 	f = 1/(2*a)*...
%<>% 	    exp(-1/2*(v/a + mu)).*...
%<>% 	    (v/(a*mu)).^(N/4 - 1/2).*...
%<>% 	    besseli(N/2 - 1,sqrt(v*mu/a));
%<>% 	
%<>% 	
%<>% 	% % then, let us compute the distribution of the inverse of the variance
%<>% 	% iv = [.0001:.00001:.01];
%<>% 	% g = 1./(iv.^2).*1/(2*a).*...
%<>% 	%     exp(-1/2*(1./(iv*a) + mu)).*...
%<>% 	%     (1./(iv*a*mu)).^(N/4 - 1/2).*...
%<>% 	%     besseli(N/2 - 1, sqrt(mu./(iv*a)));
%<>% 	% 
%<>% 	% figure(2)
%<>% 	% plot(iv,g)
%<>% 	% 
%<>% 	
%<>% 	% 
%<>% 	% % finally, the distribution of the filter coefficients
%<>% 	% b = sigma^2;
%<>% 	% w = [-10:.001:1-.01];
%<>% 	% h = (1/(2*a*b))*1./(((1-w)/b).^2).*...
%<>% 	%     exp(-1/2*(b./(a*(1-w)) + mu)).*...
%<>% 	%     (1./((1-w)/b*a*mu)).^(N/4 - 1/2).*...
%<>% 	%     besseli(N/2-1, sqrt(mu./((1-w)/b*a)));
%<>% 	% 
%<>% 	% figure(3)
%<>% 	% plot(w,h)
%<>% 	
%<>% 	
%<>% 	% and this is the equivalent expression computed from the notes
%<>% 	w = [-10:.001:1-.01];
%<>% 	h = N./(2*(1-w).^2).*...
%<>% 	    exp(-N/2*(1./(1-w) + snr)).*...
%<>% 	    (1./((1-w)*snr)).^(N/4 - 1/2).*...
%<>% 	    besseli(N/2-1, sqrt(N^2*snr./(1-w)));
%<>% 	
%<>% 	f(isnan(f)) = 0;
%<>% 	h(isnan(h)) = 0;
%<>% 	
%<>% 	% compute expected value of weight
%<>% 	Ew = sum(h.*w)*.001;
%<>% 	
%<>% 	ww = Ew;
%<>% 	hEw = N./(2*(1-ww).^2).*...
%<>% 	    exp(-N/2*(1./(1-ww) + snr)).*...
%<>% 	    (1./((1-ww)*snr)).^(N/4 - 1/2).*...
%<>% 	    besseli(N/2-1, sqrt(N^2*snr./(1-ww)));
%<>% 	
%<>% 	% compute weight corresponding to expected value of power
%<>% 	Ev  = sum(f.*v)*0.1
%<>% 	wEv = 1 - sigma^2/Ev;
%<>% 	
%<>% 	ww = wEv;
%<>% 	hwEv = N./(2*(1-ww).^2).*...
%<>% 	    exp(-N/2*(1./(1-ww) + snr)).*...
%<>% 	    (1./((1-ww)*snr)).^(N/4 - 1/2).*...
%<>% 	    besseli(N/2-1, sqrt(N^2*snr./(1-ww)));
%<>% 	
%<>% 	% gausian approximation
%<>% 	meanw = wEv;
%<>% 	varw = 2/N*(1 + 2*snr)/(1 + snr)^4;
%<>% 	
%<>% 	hg = 1/sqrt(varw*2*pi)*exp(-(w - meanw).^2/2/varw);
%<>% 	
%<>% 	% figure(1)
%<>% 	% area(v,20*f,'FaceColor','c')
%<>% 	% %alpha(.5)
%<>% 	% hold on
%<>% 	% plot(v,1-sigma^2./v,'k:')
%<>% 	% xlim([0 500])
%<>% 	% title('pdf of sample variance estimate')
%<>% 	%
%<>% 	% area(20*h,w,'FaceColor','m')
%<>% 	% %alpha(.5)
%<>% 	% plot([Ev Ev],[-1 1]*1e8,'k--')
%<>% 	% plot([-1 1]*1e8,[wEv wEv],'k--')
%<>% 	% plot([-1 1]*1e8,[Ew Ew],'k--')
%<>% 	% plot([-1 1]*1e8,[0 0],'k-','LineWidth',3)
%<>% 	% plot([0 0],[-1 1]*1e8,'k-','LineWidth',3)
%<>% 	% axis([-10 Ev + 3* sqrt(2*sigma^4/N*(1 + 2*snr)) -2 1])
%<>% 	
%<>% 	figure(3)
%<>% 	if ii == 1,
%<>% 		plot([wEv wEv],[-1 1]*1e8,'k--')
%<>% 		grid on
%<>% 		hold on
%<>% 	end
%<>% 	color = [(ii-1) length(Ns)-(ii-1) length(Ns)-(ii-1) ]/length(Ns);
%<>% 	hplots(ii) = plot(w,h,'Color',color,'LineWidth',3)
%<>% 	plot(w,hg,'--','Color',color,'LineWidth',1)
%<>% 	plot(Ew,hEw,'s','MarkerEdgeColor','k','MarkerFaceColor',color,'MarkerSize',8,'LineWidth',1)
%<>% 	maxh = max(maxh,max(h))
%<>% 	axis([-3 1 -.3 maxh+.3])
%<>% %	title('pdf of empirical wiener coefficient')
%<>% 	legends{ii} = sprintf('$N = %3d$',Ns(ii));
%<>% 	
%<>% 	
%<>% 	% % aggregation --> iterative convolution!
%<>% 	% figure(10)
%<>% 	% hold on
%<>% 	% M = 2^8;
%<>% 	% h1 = h;
%<>% 	% for i = 1:log2(M),
%<>% 	% 	h1 = conv(h1,h1,'full')*.001;
%<>% 	% 	h1 = 2*h1(1:2:end);
%<>% 	% 
%<>% 	% 	plot(w,h1)
%<>% 	% 	title(i)
%<>% 	% 	pause(.1)
%<>% 	% end
%<>% 
%<>% end
%<>% 
%<>% figure(3), hold off
%<>% hleg = legend(hplots,legends);
%<>% set(hleg,'Interpreter','latex')
%<>% legend




% parameters
snr = .05;
sigma = 20;
N = 40;

maxh = 0;
snrs = 2.^[-4 -2 0 1 2 3];
hplots = zeros(1,length(snrs));
legends = cell(1,length(snrs));
for ii = 1:length(snrs), snr = snrs(ii); 

	% first, a scaled non-central chi square
	
	v = [0:.1:8000];
	
	mu = N*snr;
	a = sigma^2/N;
	
	
	f = 1/(2*a)*...
	    exp(-1/2*(v/a + mu)).*...
	    (v/(a*mu)).^(N/4 - 1/2).*...
	    besseli(N/2 - 1,sqrt(v*mu/a));
	
	
	% % then, let us compute the distribution of the inverse of the variance
	% iv = [.0001:.00001:.01];
	% g = 1./(iv.^2).*1/(2*a).*...
	%     exp(-1/2*(1./(iv*a) + mu)).*...
	%     (1./(iv*a*mu)).^(N/4 - 1/2).*...
	%     besseli(N/2 - 1, sqrt(mu./(iv*a)));
	% 
	% figure(2)
	% plot(iv,g)
	% 
	
	% 
	% % finally, the distribution of the filter coefficients
	% b = sigma^2;
	% w = [-10:.001:1-.01];
	% h = (1/(2*a*b))*1./(((1-w)/b).^2).*...
	%     exp(-1/2*(b./(a*(1-w)) + mu)).*...
	%     (1./((1-w)/b*a*mu)).^(N/4 - 1/2).*...
	%     besseli(N/2-1, sqrt(mu./((1-w)/b*a)));
	% 
	% figure(3)
	% plot(w,h)
	
	
	% and this is the equivalent expression computed from the notes
	w = [-10:.001:1-.01];
	h = N./(2*(1-w).^2).*...
	    exp(-N/2*(1./(1-w) + snr)).*...
	    (1./((1-w)*snr)).^(N/4 - 1/2).*...
	    besseli(N/2-1, sqrt(N^2*snr./(1-w)));
	
	f(isnan(f)) = 0;
	h(isnan(h)) = 0;
	
	% compute expected value of weight
	Ew = sum(h.*w)*.001;
	
	ww = Ew;
	hEw = N./(2*(1-ww).^2).*...
	    exp(-N/2*(1./(1-ww) + snr)).*...
	    (1./((1-ww)*snr)).^(N/4 - 1/2).*...
	    besseli(N/2-1, sqrt(N^2*snr./(1-ww)));
	
	% compute weight corresponding to expected value of power
	Ev  = sum(f.*v)*0.1;
	wEv = 1 - sigma^2/Ev;
	
	ww = wEv;
	hwEv = N./(2*(1-ww).^2).*...
	    exp(-N/2*(1./(1-ww) + snr)).*...
	    (1./((1-ww)*snr)).^(N/4 - 1/2).*...
	    besseli(N/2-1, sqrt(N^2*snr./(1-ww)));
	
	% gausian approximation
	meanw = wEv;
	varw = 2/N*(1 + 2*snr)/(1 + snr)^4;
	
	hg = 1/sqrt(varw*2*pi)*exp(-(w - meanw).^2/2/varw);
	
	% figure(1)
	% area(v,20*f,'FaceColor','c')
	% %alpha(.5)
	% hold on
	% plot(v,1-sigma^2./v,'k:')
	% xlim([0 500])
	% title('pdf of sample variance estimate')
	%
	% area(20*h,w,'FaceColor','m')
	% %alpha(.5)
	% plot([Ev Ev],[-1 1]*1e8,'k--')
	% plot([-1 1]*1e8,[wEv wEv],'k--')
	% plot([-1 1]*1e8,[Ew Ew],'k--')
	% plot([-1 1]*1e8,[0 0],'k-','LineWidth',3)
	% plot([0 0],[-1 1]*1e8,'k-','LineWidth',3)
	% axis([-10 Ev + 3* sqrt(2*sigma^4/N*(1 + 2*snr)) -2 1])
	
	figure(3)
	grid on
	hold on

	color = [(ii-1) length(snrs)-(ii-1) length(snrs)-(ii-1) ]/length(snrs);
	plot([wEv wEv],[-1 1]*1e8,'--','Color',color,'LineWidth',2)
	hplots(ii) = plot(w,h,'Color',color,'LineWidth',3);
	plot(w,hg,'--','Color',color,'LineWidth',1)
	plot(Ew,hEw,'s','MarkerEdgeColor','k','MarkerFaceColor',color,'MarkerSize',8,'LineWidth',1)
	maxh = max(maxh,max(h));
	axis([-1 1 -.3 maxh+.3])
%	title('pdf of empirical wiener coefficient')
	legends{ii} = sprintf('$\\hat{\\textnormal{snr}} = %g $',snrs(ii));
	
	
	% % aggregation --> iterative convolution!
	% figure(10)
	% hold on
	% M = 2^8;
	% h1 = h;
	% for i = 1:log2(M),
	% 	h1 = conv(h1,h1,'full')*.001;
	% 	h1 = 2*h1(1:2:end);
	% 
	% 	plot(w,h1)
	% 	title(i)
	% 	pause(.1)
	% end

end

figure(3), hold off
hleg = legend(hplots,legends,'Location','Northwest');
set(hleg,'Interpreter','latex')
legend
