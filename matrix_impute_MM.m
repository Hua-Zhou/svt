function [Y,stats] = matrix_impute_MM(X,lambda,varargin)
% MATRIX_IMPUTE Impute a set of matrices using Nesterov method
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
argin.addParamValue('QN', 0, @isnumeric);
argin.addParamValue('MaxIter', 1000, @(x) isnumeric(x) && x>0);
argin.addParamValue('TolFun', 1e-5, @(x) isnumeric(x) && x>0);
argin.addParamValue('Y0', [], @(x) isnumeric(x) || isempty(x));
argin.parse(X,lambda,varargin{:});
Display = argin.Results.Display;
MaxIter = argin.Results.MaxIter;
%QN = argin.Results.QN;
TolFun = argin.Results.TolFun;
Y0 = argin.Results.Y0;

% check dimensions and retrieve mising entry information
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
if (smax<=lambda/n)
    stats.iterations = 0;
    stats.objval = norm(Wts.*Xavg,'fro')^2/2;
    stats.rank = 0;
    stats.maxlambda = smax;
    return;
end

% accumulate the first QN differences for Quasi-Newton acceleration
% if (QN>0)
%     Uqn = zeros(p1*p2,QN);
%     Vqn = zeros(p1*p2,QN);
%     for i=1:QN
%         Y_old = Y;
%         M = Xavg + W.*Y;
%         [U,s,V] = svt(M,lambda/n);
%         Y = bsxfun(@times,U,s')*V';
%         Uqn(:,i) = Y(:)-Y_old(:);
%     end
%     Vqn(:,1:QN-1) = Uqn(:,2:QN);
%     Y_old = Y;
%     M = Xavg + W.*Y;
%     [U,s,V] = svt(M,lambda/n);
%     Y = bsxfun(@times,U,s')*V';
%     Vqn(:,QN) = Y(:)-Y_old(:);
%     old_secant = 1;
%     C = Uqn'*(Uqn-Vqn);
% end

% main loop
objval = inf;
for iter=1:MaxIter

    % thesholding intermediate matrix
    %M = Xavg + W.*Y;
    WY = Wts.*Y;  
    D = Xavg-WY;
    [U,s,V] = svt(@MAtimesVec,'m',p1,'n',p2,'lambda',lambda/n); % call svt
    %[U,s,V] = svt(M,'lambda',lambda/n); % call svt
    %[U,s,V] = fsvt(M,lambda/n);                 % call fsvt
    s = diag(s)-lambda/n;                              % shrinkage
    %Y_old = Y;
    Y = bsxfun(@times,U,s')*V';
    objval_old = objval;
    
    %if (QN==0)      % no Quasi-Newton acceleration
    objval = norm(Wts.*(Y-Xavg),'fro')^2/2 + lambda*sum(s);
%     elseif (QN>0)   % Quasi-Newton acceleration
%         % do one more MM step to accumulate secant pairs
%         Uqn(:,old_secant) = Y(:) - Y_old(:);
%         M = Xavg + W.*Y;
%         [U,s,V] = svt(M,lambda/n);
%         Y_old = Y;
%         Y = bsxfun(@times,U,s')*V';
%         Vqn(:,old_secant) = Y(:) - Y_old(:);
%         C(old_secant,:) = Uqn(:,old_secant)'*(Uqn-Vqn);
%         C(:,old_secant) = Uqn'*(Uqn(:,old_secant)-Vqn(:,old_secant));
%         new_secant = old_secant;
%         old_secant = mod(old_secant,QN)+1;
%         objval_MM = norm(Wts.*(Y-Xavg),'fro')^2/2 + lambda*sum(s);
%         % quasi-Newton jump
%         Y_QN = Y_old + reshape(Vqn*(C\(Uqn'*Uqn(:,new_secant))),p1,p2);
%         M = Xavg + W.*Y_QN;
%         [U,s,V] = svt(M,lambda/n);
%         Y_QN = bsxfun(@times,U,s')*V';
%         objval_QN = norm(Wts.*(Y_QN-Xavg),'fro')^2/2 + lambda*sum(s);
%         % choose MM vs QN jump
%         if (objval_QN<objval_MM)
% %             disp('yeah!');
%             Y = Y_QN;
%             objval = objval_QN;
%         else
%             objval = objval_MM;
%         end
%     end
    
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
function MAvec = MAtimesVec(vec, varargin)
    argin = inputParser;
    argin.addRequired('vec');
    argin.addOptional('trans', false, @islogical);
    argin.parse(vec,varargin{:});

    trans = argin.Results.trans;

    if trans
       MAvec = (vec'*D)' + (vec'*Y)';
    else
       MAvec = D*vec + Y*vec;
    end
    
end

end