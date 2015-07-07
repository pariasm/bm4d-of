% data


prodmatrix2_eigs  = 0.0285;
prodmatrix2_idsvd = 0.1392;

matrixsvd2_idsvd = 0.6220; % * step2_lr
matrixeig2_eigs  = 0.6338; % * step2_fr

%
% (time matrix eigs) / (time matrix inv)
%
covdecomp2_idsvd_over_eigs = prodmatrix2_eigs / prodmatrix2_idsvd * matrixsvd2_idsvd / matrixeig2_eigs;
disp(['s10 | Reduction in computation of principal directions:' ])
disp(covdecomp2_idsvd_over_eigs)

%
% (time step 2 low rank) / (time step 2 full rank)
%
step2_lr_over_fr = prodmatrix2_eigs / prodmatrix2_idsvd;
disp(['Overall reduction in step 2 bayes estimate:' ])
disp(step2_lr_over_fr)


