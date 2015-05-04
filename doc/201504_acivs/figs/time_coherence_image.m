v1 = mean(single(imread('/data/pablo/boucantrin/denoising/projects/video_nlbayes/results/vnlbayes/table_1ps5_2ps5_2wx37_2np93/tennis_s40_tr2/deno_121.png')),3);
v2 = mean(single(imread('/data/pablo/boucantrin/denoising/projects/video_nlbayes/results/vnlbayes/table_1ps5_2ps5_2wx37_2np93/tennis_s40_tr2/deno_122.png')),3);
v3 = mean(single(imread('/data/pablo/boucantrin/denoising/projects/video_nlbayes/results/vnlbayes/table_1ps5_2ps5_2wx37_2np93/tennis_s40_tr2/deno_123.png')),3);
v4 = mean(single(imread('/data/pablo/boucantrin/denoising/projects/video_nlbayes/results/vnlbayes/table_1ps5_2ps5_2wx37_2np93/tennis_s40_tr2/deno_124.png')),3);
v5 = mean(single(imread('/data/pablo/boucantrin/denoising/projects/video_nlbayes/results/vnlbayes/table_1ps5_2ps5_2wx37_2np93/tennis_s40_tr2/deno_125.png')),3);

td1 = abs(v1 - v2) + ...
      abs(v2 - v2) + ...
      abs(v3 - v3) + ...
      abs(v4 - v5) / 4;

td1 = min(td1/20*255,255);

u1 = mean(single(imread('/data/pablo/purple/dev/results/VBM4D/tennis_s40/deno_121.png')),3);
u2 = mean(single(imread('/data/pablo/purple/dev/results/VBM4D/tennis_s40/deno_122.png')),3);
u3 = mean(single(imread('/data/pablo/purple/dev/results/VBM4D/tennis_s40/deno_123.png')),3);
u4 = mean(single(imread('/data/pablo/purple/dev/results/VBM4D/tennis_s40/deno_124.png')),3);
u5 = mean(single(imread('/data/pablo/purple/dev/results/VBM4D/tennis_s40/deno_125.png')),3);

td2 = abs(u1 - u2) + ...
      abs(u2 - u2) + ...
      abs(u3 - u3) + ...
      abs(u4 - u5) / 4;

td2 = min(td2/20*255,255);

imwrite(uint8(td1),'time_diff_vnlb_s40.png');
imwrite(uint8(td2),'time_diff_bm4d_s40.png');
