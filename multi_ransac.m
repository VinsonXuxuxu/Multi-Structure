
function [M, sh_old, trialcount] = multi_ransac(x, fittingfn, distfn, degenfn, s, t, numMod, feedback, ...
                               maxDataTrials, maxTrials)

    % Test number of parameters
    error ( nargchk ( 7, 10, nargin ) );
    
    if nargin < 10; maxTrials = 1000;    end;
    if nargin < 9; maxDataTrials = 100; end;
    if nargin < 8; feedback = 0;        end;
    
    [rows, npts] = size(x);
    
    p = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers (probably should be a parameter)

    bestM = NaN;      % Sentinel value allowing detection of solution failure.
    trialcount = 0;
    bestscore =  0;
    N = 1;            % Dummy initialisation for number of trials.
    
    sh_old = cell(1, numMod);
    sz_sh_old = zeros(1, numMod); 
    while N > trialcount
        
        % Select at random s datapoints to form a trial model, M.
        % In selecting these points we have to check that they are not in
        % a degenerate configuration.
        idx = 1:npts; 
        sh = cell(1, numMod); 
        sz_sh = zeros(1, numMod); 
        
        for i = 1:numMod
            degenerate = 1;
            count = 1;
            while degenerate
                % Generate s random indicies in the range 1..npts
                % (If you do not have the statistics toolbox with randsample(),
                % use the function RANDOMSAMPLE from my webpage)
    %             if ~exist('randsample', 'file')
                    ind = idx(randomsample(numel(idx), s));

                % Test that these points are not a degenerate configuration.
                degenerate = feval(degenfn, x(:,ind));

                if ~degenerate
                    % Fit model to this random selection of data points.
                    % Note that M may represent a set of models that fit the data in
                    % this case M will be a cell array of models
                    M = feval(fittingfn, x(:,ind));

                    % Depending on your problem it might be that the only way you
                    % can determine whether a data set is degenerate or not is to
                    % try to fit a model and see if it succeeds.  If it fails we
                    % reset degenerate to true.
                    if isempty(M)
                        degenerate = 1;
                    end
                end

                % Safeguard against being stuck in this loop forever
                count = count + 1;
                if count > maxDataTrials
                    warning('Unable to select a nondegenerate data set');
                    break
                end
            end

            [inliers, M] = feval(distfn, M, x(:, idx), t);
            sh{i} = idx(inliers); 
            idx(inliers) = []; 
            % Find the number of inliers to this model.
            ninliers = length(inliers);
            sz_sh(i) = ninliers; 
        end
        
        sh_all = [sh_old, sh]; 
        sz_all = [sz_sh_old, sz_sh]; 
        sh_old = cell(1, numMod); 
        sz_sh_old = zeros(1, numMod); 
        [mx, id] = sort(sz_all, 2, 'descend'); 
        count = 1; 
        cnt = 0; 
        while(count <= numMod && cnt < numel(sh_all))
            cnt = cnt + 1; 
            flag = 0; 
            for i=1:count
                if numel(intersect(sh_all{id(cnt)}, sh_old{i})) > 0; 
                    flag = 1; 
                end
            end
            if ~flag
                sh_old{count} = sh_all{id(cnt)}; 
                sz_sh_old(count) = numel(sh_old{count}); 
                count = count + 1; 
            end
        end
            
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  min(sz_sh_old)/npts;
            pNoOutliers = 1 -  fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
        
        trialcount = trialcount+1;
        if feedback
            fprintf('trial %d out of %d         \r',trialcount, ceil(N));
        end

        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            warning( ...
            sprintf('ransac reached the maximum number of %d trials',...
                    maxTrials));
            break
        end
    end
    
    M = zeros(1, npts); 
    
    for i=1:numMod
        M(sh_old{i}) = i;
    end
    
end
