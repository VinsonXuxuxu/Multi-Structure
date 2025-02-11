
%% Multi-Structure Fitting Test On Star-Shape Data

% This code just simply run the T-Linkage algorithm on the example data set
% "star5".
% Loading data: X contains data points, whereas G is the ground truth
% segmentation
load 'star5.mat'; N = size(X,2);
% Model specification: T-Linkage need to define a distance (distFun)
% between points and models.
% hpFun given cardmss points returns an estimate model, fit_model is
% a function for least square fitting.
% In this example we want to estimate lines so distFun is the euclidean
% distance between a point from a line in the plane and cardmss=2.
% Other  possible models are 'line', 'circle',
% fundamental matrices ('fundamental') and 'subspace4' (look in 'model_spec' folder).
%

numMod = 5; 

[distFun, hpFun, fit_model, cardmss] = set_model('line');
%% Conceptual representation of points

%T-linkage starts, as Ransac with random sampling:
% Unform sampling can be adopted
S = mssUniform(N, 5*N,cardmss); 
% in order to reduce the number of hypotheses also a localized sampling can
% be used:
%
%       D = pdist(X','euclidean');  D = squareform(D);
%       S = mssNorm( X, D, 2*N, cardmss);
%

H = hpFun( X, S ); 
% generating a pool of putative hypotheses H.
% The residuals R between points and model
R = res( X, H, distFun );
% are used for representing points in a conceptual space.
% In particular a preference matrix P is built depicting by rows points
% preferences.
% 

epsilon= 1.0e-1; %An inlier threshold value  epsilon has to be specified.

CT = T_linkage(R, epsilon, cardmss); 
MT = fit_model( X, CT); 
clrsT = distinguishable_colors(size(MT, 2)+1); 

CJ = J_linkage(R, epsilon, cardmss); 
MJ = fit_model( X, CJ); 
clrsJ = distinguishable_colors(size(MJ, 2)+1); 

CM = M_shift(X, H, R, epsilon, distFun); 
MM = fit_model( X, CM); 
clrsM = distinguishable_colors(size(MM, 2)+1); 

CRS = Ransac_seq(X, epsilon, numMod); 
MRS = fit_model( X, CRS); 
clrsRS = distinguishable_colors(size(MRS, 2)+1); 

CRM = Ransac_multi(X, epsilon, numMod); 
MRM = fit_model( X, CRM); 
clrsRM = distinguishable_colors(size(MRM, 2)+1); 

%% Showing results

MG = fit_model( X, G); 
clrsG = distinguishable_colors(size(MG, 2)+1); 

figure
subplot(2,3,1); gscatter(X(1,:),X(2,:),C, clrsG); lineplot(MG,clrsG);  axis equal; title('GroundTruth'); legend off
subplot(2,3,2); gscatter(X(1,:),X(2,:),CT, clrsT); lineplot(MT,clrsT); axis equal; title('T-linkage'); legend off
subplot(2,3,3); gscatter(X(1,:),X(2,:),CJ, clrsJ); lineplot(MJ,clrsJ); axis equal; title('J-linkage'); legend off
subplot(2,3,4); gscatter(X(1,:),X(2,:),CM, clrsM); lineplot(MM,clrsM); axis equal; title('Mean-Shift'); legend off
subplot(2,3,5); gscatter(X(1,:),X(2,:),CRS, clrsRS); lineplot(MRS,clrsRS); axis equal; title('Sequential-Ransac'); legend off
subplot(2,3,6); gscatter(X(1,:),X(2,:),CRM, clrsRM); lineplot(MRM,clrsRM); axis equal; title('Multi-Ransac'); legend off

%% Reference
% When using the code in your research work, please cite the following paper:
% Luca Magri, Andrea Fusiello, T-Linkage: A Continuous Relaxation of
% J-Linkage for Multi-Model Fitting, CVPR, 2014.
%
% For any comments, questions or suggestions about the code please contact
% luca (dot) magri (at) unimi (dot) it



