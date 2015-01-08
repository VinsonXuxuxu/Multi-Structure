
function C = Ransac_seq(X, epsilon, numMod)

    N = size(X, 2); 
    X = [X; ones(1, N)]; 
    C = zeros(N, 1); 

    feedback = 0; 
    fit_model = @myFitLine;
    distFun    = @lineptdist;
    degenfn   = @isdegenerate;
    cardmss = 2; 
    
    idx = 1:N; 
    for i=1:numMod
        [L, inliers, trialcount] = ransac(X(:, idx), fit_model, distFun, degenfn, cardmss, epsilon, feedback);
        C(idx(inliers)) = i; 
        idx(inliers) = []; 
    end

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