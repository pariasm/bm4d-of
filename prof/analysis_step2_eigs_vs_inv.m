% data


% s10 -------------------------------------------------------------------
cov2_fr = 0.2094;
cov2_lr = 0.3409;

covinv2_fr = 0.4139; % * step2_fr
coveig2_lr = 0.6220; % * step2_lr

mprod2_fr = 1 - cov2_fr - covinv2_fr;
mprod2_lr = 1 - cov2_lr - coveig2_lr;

%
% (time matrix eigs) / (time matrix inv)
%
covdecomp2_lr_over_fr = cov2_fr / cov2_lr * coveig2_lr / covinv2_fr;
disp(['s40 | Reduction in cov matrix decomposition:' ])
disp(covdecomp2_lr_over_fr)

%
% (time matrix prod low rank) / (time matrix prod full rank)
%
matrixprod2_lr_over_fr = cov2_fr / cov2_lr * mprod2_lr / mprod2_fr;
disp(['s40 | Reduction in matrix product for group denoising:' ])
disp(matrixprod2_lr_over_fr)

%
% (time step 2 low rank) / (time step 2 full rank)
%
step2_lr_over_fr = cov2_fr / cov2_lr;
disp(['s40 | Overall reduction in step 2 bayes estimate:' ])
disp(step2_lr_over_fr)




% s20 -------------------------------------------------------------------
cov2_fr = 0.2094;
cov2_lr = 0.3409;

covinv2_fr = 0.4139; % * step2_fr
coveig2_lr = 0.6220; % * step2_lr

mprod2_fr = 1 - cov2_fr - covinv2_fr;
mprod2_lr = 1 - cov2_lr - coveig2_lr;

%
% (time matrix eigs) / (time matrix inv)
%
covdecomp2_lr_over_fr = cov2_fr / cov2_lr * coveig2_lr / covinv2_fr;
disp(['s40 | Reduction in cov matrix decomposition:' ])
disp(covdecomp2_lr_over_fr)

%
% (time matrix prod low rank) / (time matrix prod full rank)
%
matrixprod2_lr_over_fr = cov2_fr / cov2_lr * mprod2_lr / mprod2_fr;
disp(['s40 | Reduction in matrix product for group denoising:' ])
disp(matrixprod2_lr_over_fr)

%
% (time step 2 low rank) / (time step 2 full rank)
%
step2_lr_over_fr = cov2_fr / cov2_lr;
disp(['s40 | Overall reduction in step 2 bayes estimate:' ])
disp(step2_lr_over_fr)




% s40 -------------------------------------------------------------------
cov2_fr = 0.2094;
cov2_lr = 0.3409;

covinv2_fr = 0.4139; % * step2_fr
coveig2_lr = 0.6220; % * step2_lr

mprod2_fr = 1 - cov2_fr - covinv2_fr;
mprod2_lr = 1 - cov2_lr - coveig2_lr;

%
% (time matrix eigs) / (time matrix inv)
%
covdecomp2_lr_over_fr = cov2_fr / cov2_lr * coveig2_lr / covinv2_fr;
disp(['s40 | Reduction in cov matrix decomposition:' ])
disp(covdecomp2_lr_over_fr)

%
% (time matrix prod low rank) / (time matrix prod full rank)
%
matrixprod2_lr_over_fr = cov2_fr / cov2_lr * mprod2_lr / mprod2_fr;
disp(['s40 | Reduction in matrix product for group denoising:' ])
disp(matrixprod2_lr_over_fr)

%
% (time step 2 low rank) / (time step 2 full rank)
%
step2_lr_over_fr = cov2_fr / cov2_lr;
disp(['s40 | Overall reduction in step 2 bayes estimate:' ])
disp(step2_lr_over_fr)

