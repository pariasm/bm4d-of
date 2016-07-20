% loads a video and allows the user to select a pixel
% for the selected pixel, it displays the locations 
% of the most similar patches

%clear all


% parameters
prms1.wx = 25;
prms1.wt = 0;
prms1.px = 5;
prms1.pt = 3;
prms1.np = 10;
prms1.r  = -1;

prms2.wx = 25;
prms2.wt = 2;
prms2.px = 5;
prms2.pt = 3;
prms2.np = 10;
prms2.r  = -1;

sigma = 40;

% binary that computes patch distances
pgbin = '/home/pariasm/Work/denoising/projects/video_nlbayes3d/build/bin/patch_group';

% path to data
base_path = '../../../';
orig_path = [base_path 'data/derf/'];

% sequence data
seqnames = {'bus','coastguard','foreman','tennis'};
first    = {1    , 1          , 1       , 1      };
last     = {150  , 300        , 300     , 150    };

sequences = struct('name',seqnames,'first',first,'last',last);
clear seqnames first last

% select a sequence
seq = sequences(1);

ori_pat = [orig_path seq.name '/%03d.png'];                     % ground truth

% load
orig = single(load_sequence(ori_pat, seq.first, seq.last));
%randn('seed',0)
%nisy = orig + sigma*randn(size(orig));

% process only luminance
%nisy = mean(nisy,3);
%orig = 128 + 0*mean(orig,3);
orig = mean(orig,3);
%sigma = sigma/sqrt(3);

%randn('seed',0)
nisy = orig + sigma*randn(size(orig));

% crop
cropx = [170,270];
cropy = [ 30,130];
cropt = [  1, 10];

nisy = nisy(cropy(1):cropy(2),cropx(1):cropx(2),:,cropt(1):cropt(2));
orig = orig(cropy(1):cropy(2),cropx(1):cropx(2),:,cropt(1):cropt(2));


%prms1.filter_type = 'pos';
%[deno1_pos, aggw1, bias1_pos] = nlbayes_step(nisy, []       , orig, sigma, prms1);
%[deno2_pos, aggw2, bias2_pos] = nlbayes_step(nisy, deno1_pos, orig, sigma, prms2);
%
%o = orig(5:end-5,5:end-5,:,:);
%b1p = bias1_pos(5:end-5,5:end-5,:,:);
%b2p = bias2_pos(5:end-5,5:end-5,:,:);
%d1p = deno1_pos(5:end-5,5:end-5,:,:);
%d2p = deno2_pos(5:end-5,5:end-5,:,:);
%disp(20*log10(255/sqrt(norm(o(:) - d1p(:))^2/length(o(:)))))
%disp(20*log10(255/sqrt(norm(o(:) - d2p(:))^2/length(o(:)))))

%[deno3_pos, aggw3, bias3_pos] = nlbayes_step(nisy, deno2_pos, orig, sigma, prms2);
%[deno4_pos, aggw4, bias4_pos] = nlbayes_step(nisy, deno3_pos, orig, sigma, prms2);
%[deno5_pos, aggw5, bias5_pos] = nlbayes_step(nisy, deno4_pos, orig, sigma, prms2);
%d3p = deno3_pos(5:end-5,5:end-5,:,:);
%d4p = deno4_pos(5:end-5,5:end-5,:,:);
%d5p = deno5_pos(5:end-5,5:end-5,:,:);
%disp(20*log10(255/sqrt(norm(o(:) - d3p(:))^2/length(o(:)))))
%disp(20*log10(255/sqrt(norm(o(:) - d4p(:))^2/length(o(:)))))
%disp(20*log10(255/sqrt(norm(o(:) - d5p(:))^2/length(o(:)))))

prms1.filter_type = 'neg';
%[deno1_neg, aggw1, bias1_neg] = nlbayes_step(nisy, []       , orig, sigma, prms1);
%[deno2_neg, aggw2, bias2_neg] = nlbayes_step(nisy, deno1_neg, orig, sigma, prms2);
[deno1_neg, aggw1, bias1_neg] = nlbayes_step(nisy, []       , [], sigma, prms1);
[deno2_neg, aggw2, bias2_neg] = nlbayes_step(nisy, deno1_neg, [], sigma, prms2);

b1n = bias1_neg(5:end-5,5:end-5,:,:);
b2n = bias2_neg(5:end-5,5:end-5,:,:);
d1n = deno1_neg(5:end-5,5:end-5,:,:);
d2n = deno2_neg(5:end-5,5:end-5,:,:);
disp(20*log10(255/sqrt(norm(o(:) - d1n(:))^2/length(o(:)))))
disp(20*log10(255/sqrt(norm(o(:) - d2n(:))^2/length(o(:)))))

%[deno3_neg, aggw3, bias3_neg] = nlbayes_step(nisy, deno2_neg, orig, sigma, prms2);
%[deno4_neg, aggw4, bias4_neg] = nlbayes_step(nisy, deno3_neg, orig, sigma, prms2);
%[deno5_neg, aggw5, bias5_neg] = nlbayes_step(nisy, deno4_neg, orig, sigma, prms2);
%d3n = deno3_neg(5:end-5,5:end-5,:,:);
%d4n = deno4_neg(5:end-5,5:end-5,:,:);
%d5n = deno5_neg(5:end-5,5:end-5,:,:);
%disp(20*log10(255/sqrt(norm(o(:) - d3n(:))^2/length(o(:)))))
%disp(20*log10(255/sqrt(norm(o(:) - d4n(:))^2/length(o(:)))))
%disp(20*log10(255/sqrt(norm(o(:) - d5n(:))^2/length(o(:)))))


%prms1.filter_type = 'neg-inv';
%[deno1_ngi, aggw1] = nlbayes_step(nisy, []       , orig, sigma, prms1);
%[deno2_ngi, aggw2] = nlbayes_step(nisy, deno1_ngi, orig, sigma, prms2);

%% romano's iteration
%prms2.filter_type = 'pos';
%[deno1_rom, aggw1] = nlbayes_step(0.5*(nisy + deno2_neg), []       , orig, 0.5*sigma, prms1);
%[deno2_rom, aggw2] = nlbayes_step(0.5*(nisy + deno2_neg), deno1_rom, orig, 0.5*sigma, prms2);
%deno2_rom = deno2_rom - 0.5*deno2_neg;

%d1i = deno1_ngi(5:end-5,5:end-5,:,:);
%d2i = deno2_ngi(5:end-5,5:end-5,:,:);

%disp(20*log10(255/sqrt(norm(o(:) - d1i(:))^2/length(o(:)))))
%disp(20*log10(255/sqrt(norm(o(:) - d2i(:))^2/length(o(:)))))

