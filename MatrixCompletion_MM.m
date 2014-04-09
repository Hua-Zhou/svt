function [Y,stats] = MatrixCompletion_MM(X,lambda,varargin)
% MATRIX_IMPUTE Impute a set of matrices using MM algorithm
%
% INPUT:
%   X - p1-by-p2-by-n data matrices for imputation; missing values are nan
%
% Output:
%   Y - the imputation matrix
%   stats - algorithmic statistics
%
% COPYRIGHT North Carolina State University 
% AUTHOR: Hua Zhou (hua_zhou@ncsu.edu)

% parse inputs
argin = inputParser;
argin.addRequired('X', @isnumeric);
argin.addRequired('lambda', @isnumeric);
argin.addParamValue('Display', 'iter', @(x) ischar(x));
argin.addParamValue('MaxIter', 1000, @(x) isnumeric(x) && x>0);
argin.addParamValue('TolFun', 1e-5, @(x) isnumeric(x) && x>0);
argin.addParamValue('Y0', [], @(x) isnumeric(x) || isempty(x));
argin.addParamValue('method','stru_svt',@(x) ischar(x));
argin.parse(X,lambda,varargin{:});
Display = argin.Results.Display;
MaxIter = argin.Results.MaxIter;
TolFun = argin.Results.TolFun;
Y0 = argin.Results.Y0;
method = argin.Results.method;

% check dimensions and retrieve missing entry information
[p1,p2,n] = size(X);
Wts = n-sum(isnan(X),3);    % weight matrix in objective
W = 1-sum(~isnan(X),3)/n;   % weight matrix for MM algorithm
Xavg = mean(X,3);
Xavg(isnan(Xavg)) = 0;

% initialize
if (isempty(Y0))
    Y = zeros(p1,p2);
else
    Y = Y0;
end

% output max lambda if get all 0 singular values
smax = svds(Xavg + W.*Y,1);
if (smax-lambda/n<=1e-5)               %Tol of difference
    stats.iterations = 0;
    stats.objval = norm(Wts.*Xavg,'fro')^2/2;
    stats.rank = 0;
    stats.maxlambda = smax;
    return;
end

% main loop
objval = inf;
for iter=1:MaxIter
    % thesholding intermediate matrix
    if (strcmpi(method,'stru_svt'))
        D = Xavg-Wts.*Y;
        D = sparse(D);
        [U,s,V] = svt(@MAtimesVec,'m',p1,'n',p2,'lambda',lambda/n,'method','succession');
%         if iter==1
%             [U,s,V] = svt(@MAtimesVec,'m',p1,'n',p2,'lambda',lambda/n,'method','succession'); 
%         else
%             [U,s,V] = svt(@MAtimesVec,'m',p1,'n',p2,'lambda',lambda/n,...
%                 'k',lens,'method','succession');
%         end
        s = diag(s)-lambda/n;    % shrinkage
%          lens = length(s);
    elseif (strcmpi(method,'svt'))
        M = Xavg + W.*Y;
        [U,s,V] = svt(M,'lambda',lambda/n,'method','succession'); % call svt
%         if iter==1
%             [U,s,V] = svt(M,'lambda',lambda/n,'method','succession'); % call svt
%         else
%             [U,s,V] = svt(M,'lambda',lambda/n,'k',lens,'method','succession'); % call svt
%         end
        s = diag(s)-lambda/n;    % shrinkage
%         lens = length(s);
    elseif (strcmpi(method,'full'))
        M = Xavg + W.*Y;
        [U,s,V] = fsvt(M,lambda/n);                 % call fsvt
        s = s-lambda/n;    % shrinkage        
    end
    Y = bsxfun(@times,U,s')*V';
    objval_old = objval;
    objval = norm(Wts.*(Y-Xavg),'fro')^2/2 + lambda*sum(s);

    % stopping rule
    if (abs(objval_old-objval)<TolFun*(abs(objval_old)+1))
        break;
    end
    
    % display
    if (~strcmpi(Display,'off'))
        display(['iter ' num2str(iter) ', objval=' num2str(objval)]);
    end

end

% collect algorithmic statistics
stats.iterations = iter;
stats.objval = objval;
stats.rank = length(s);

% Subfunction for utilizing matrix structure of sparse plus low rank
function MAvec = MAtimesVec(vec, trans)
%     if trans
%          MAvec = (vec'*D)' + (vec'*Y)';
%     else
%          MAvec = D*vec + Y*vec;
%     end 
    if iter == 1
       if trans
         MAvec = (vec'*D)' + (vec'*Y)';
       else
         MAvec = D*vec + Y*vec;
       end
    else
       if trans
         MAvec = (vec'*D)' + V*(s.*(vec'*U)');
       else
         MAvec = D*vec + U*(s.*(V'*vec));
       end
    end
    
end

end