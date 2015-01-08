
function C = Ransac_multi(X, epsilon, numMod)

    N = size(X, 2); 
    X = [X; ones(1, N)]; 
    C = zeros(N, 1); 

    feedback = 1; 
    fit_model = @myFitLine;
    distFun    = @lineptdist;
    degenfn   = @isdegenerate;
    cardmss = 2; 
    
    [C, inliers, trialcount] = multi_ransac(X, fit_model, distFun, degenfn, cardmss, epsilon, numMod, feedback);
  
end


function [inliers, xn] = lineptdist(xn, XY, t)

    X = XY(1:end-1, :)';  
    Y = XY(end, :)'; 
    d = (X*xn - Y)/norm(xn); 
    
    inliers = find(abs(d) < t);
%     inliers = find(abs(d/sqrt(1+sum(xn(1:end-1).^2))) < t);
    
end

function r = isdegenerate(X)
    r = 0; 
%     r = norm(X(:,1) - X(:,2)) < eps;
end