%% The demo for constructing random forest based affinity matrices
% 
% @Author: Xiatian (Eddy) Zhu
% @Date: 17 June. 2014


addpath('random forest');
addpath('SPClust');

%% Load data
load('data');

%% Train a clustering random forest
% Parameters
ntree = 200;
mtry = -1;
extra_options.proximity = 1;
extra_options.nodesize = 1;

RF_model = classRF_train(X, [], ntree, mtry, extra_options);


%% Build affinity matrix
% The affinity matrix by the binary affinity model
A_Bi = RF_model.proximity;
figure(1);
imagesc(A_Bi);
title('ClustRF-Bi');

% The affinity matrix by the uniform ClustRF-Strct model
disp('To construct affinity by ClustRF-Strct(Unfm)');
A_Unfm = build_ClustRF_Strct_A(X, RF_model, 'Uniform');
figure(2);
imagesc(A_Unfm);
title('ClustRF-Strct (Unfm)');

% The affinity matrix by the adaptive ClustRF-Strct model
disp('To construct affinity by ClustRF-Strct(Adpt)');
A_Adpt = build_ClustRF_Strct_A(X, RF_model, 'Adaptive');
figure(3);
imagesc(A_Adpt);
title('ClustRF-Strct (Adpt)');


%% Perform spectral clustering
num_clst = 6;
Cl_Bi = SPClustering(A_Bi, num_clst);
Cl_Unfm = SPClustering(A_Unfm, num_clst);
Cl_Adpt = SPClustering(A_Adpt, num_clst);


%% Compare the clustering results
ARI_Bi = adjust_rand_index(Cl_Bi, Y);
ARI_Unfm = adjust_rand_index(Cl_Unfm, Y);
ARI_Adpt = adjust_rand_index(Cl_Adpt, Y);


fprintf('ARI score comparison: \nClustRF_Bi: %f \nClustRF_Strct(Unfm): %f \nClustRF_Strct(Adpt): %f\n', ARI_Bi, ARI_Unfm, ARI_Adpt);

