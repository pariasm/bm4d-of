% generate laplace distributed values
for sg = 2.^[0:.5:8],

n = 40;
mu = 10;
%sg = 5;
sg_noise = 10;

sg


a_mle_mu    = zeros(5000,1);
a_mle_sg    = zeros(5000,1);
a_mean_mu   = zeros(5000,1);
a_mean_sg   = zeros(5000,1);
b_mean_mu   = zeros(5000,1);
b_mean_sg   = zeros(5000,1);
b_median_mu = zeros(5000,1);
b_median_sg = zeros(5000,1);
for xx = 1:5000;
	alpha = laprnd(n,1,mu,sg);
	beta  = alpha + sg_noise * randn(n,1);


	% MLE estimators from clean data
	a_mle_mu(xx) = median(alpha);
	a_mle_sg(xx) = sqrt(2)*mean(abs(alpha - a_mle_mu(xx)));

	% mean and variance from clean data
	a_mean_mu(xx) = mean(alpha);
	a_mean_sg(xx) = sqrt(mean((alpha - a_mean_mu(xx)).^2));

	% estimators from noisy data
	b_mean_mu(xx) = mean(beta);
	b_mean_sg(xx) = sqrt(mean((beta - b_mean_mu(xx)).^2) - sg_noise*sg_noise);

	b_median_mu(xx) = median(beta);
	b_median_sg(xx) = sqrt(mean((beta - b_median_mu(xx)).^2) - sg_noise*sg_noise);

	% this is not actually an mle estimator, but we called mle by analogy
	% with the mle estimator from alpha
	b_mle_mu(xx) = median(beta);
	b_mle_sg(xx) = sqrt(2*mean(abs(beta - b_mle_mu(xx)))^2 - sg_noise*sg_noise);
end



%figure(1), clf, hold on
%plot(sort(abs(mu - a_mle_mu)),'k')
%plot(sort(abs(mu - a_mean_mu)),'r')
%plot(sort(abs(mu - b_median_mu)),'k--')
%plot(sort(abs(mu - b_mean_mu)),'r--')
%hold off
%axis([0 5000 0 50])
%
figure(2), clf, hold on
plot(sort(abs(sg - a_mle_sg)),'k')
plot(sort(abs(sg - a_mean_sg)),'r')
plot(sort(abs(sg - b_median_sg)),'k--')
plot(sort(abs(sg - b_mean_sg)),'r--')
plot(sort(abs(sg - b_mle_sg)),'g--')
hold off
axis([0 5000 0 50])
pause(.01)

end


