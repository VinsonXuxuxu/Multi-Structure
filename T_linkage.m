
function C = T_linkage(R, epsilon, cardmss)
P  = prefMat( R, epsilon, 1 );

%% Clustering

%T-Linkage clustering follow a bottom up scheme in the preference space

C = tlnk(P);

% C is a vector of labels, points belonging to the same models share the
% same label.
%% Outlier rejection step

%T-Linkage fit a model to all the data points. Outlier can be found in
%different ways (T-Linkage is agonostic about the outlier rejection strategy),
%for example discarding too small cluster, or exploiting the randomness of
%a model.
 C  = outlier_rejection_card( C, cardmss );
% Outliers are labelled by '0'
