
function [C, M] = PEaRL_multi(X, H, R)

    N = size(X,2);
    L = size(H, 2); 

    lambda = 5; 
    pD = (exp(-pdist2(X', X').^2/lambda^2) - eye(N)); 

    %% Clustering
    h = GCO_Create(N, L); 
    
    GCO_SetDataCost(h, R');
    GCO_SetSmoothCost(h, ones(L)); 
    GCO_SetNeighbors(h, pD); 
    GCO_SetLabelCost(h, 10); 
    
%     GCO_Expansion(h); 
    GCO_ExpandOnAlpha(h, 1:L); 
    
    C = GCO_GetLabeling(h); 
%     C  = outlier_rejection_card( C, 2 );
    M = H(:); 
    
end
