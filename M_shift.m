
function C = M_shift(X, H, R, epsilon, distFun)

N = size(X,2);
P  = prefMat( R, epsilon/2, 0 ); 
H = H(:, sum(P) > N/10); 

%% Clustering

%T-Linkage clustering follow a bottom up scheme in the preference space

bandwidth = 1; 
H = H./(ones(3, 1)*H(3, :)); 
tic
[clustCent,point2cluster, ~] = MeanShiftCluster(H, bandwidth);
toc
[~, Gl]  = outlier_rejection_card( point2cluster', 1 ); 
R = res( X, clustCent(:, Gl), distFun ); 
[mn, C] = min(R, [], 2); 

C(mn > epsilon) = 0; 

