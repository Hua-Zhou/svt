function[U,S,V,flag] = svt(A,varargin)
% svt Singular value thresholding
%   svt computes the singular values exceeding user defined threshold and
%   associated singular vectors. It can also be used for top singular value
%   decomposition. It handles sparse matrix and other sturctued matrix. In
%   the later case, user inputs a function handle instead of the data
%   matrix to exploit the matrix structure for fast matrix-vector
%   multiplication.
% 
% USAGE:
%   [U,S,V,flag] = svt(A,'PARAM1',val1,'PARAM2',val2...)
%   S = svt(A,'PARAM1',val1,'PARAM2',val2...)
%
% INPUT:
%   A - m-by-n matrix or a function handle provided by user
%   A should be large and sparse if A is a matrix
%
% OUTPUT:
%   U - left singular vectors
%   S - diagonal matrix of thresholded singular values
%   V - right singular vectors
%   flag - if 0, iterative eigs converged; 1, eigs not converged
%
% Optional parameter name/value pairs are:      
%   'lambda' - threshold (default: NaN). When the value is NaN, svt
%       retrieves top singular values/vectors.
%   'k' - number of singular values to try (default: 6)
%   'incre' - increment to catch the threshold (default: 5)
%   'm' - number of rows. Required for function handle input.
%   'n' - number of columns. Required for function handle input.
%   'tol' - convergence tolerance (default: 1e-10)
%   'maxit' - maximum number of eigs successions (default: 300) 
%   'method' - deflation method (default: 'deflation'). If specified as 
%       'succession', a succession method is applied for thresholding.
%
% Examples:
%  Singular value decomposition for the first 6 singular values. 
%  [U,S,V] = svt(A)
%
%  Singular value decomposition for the first 15 singular values.
%  [U,S,V] = svt(A,'k',15) 
%
%  Singular value thresholding, only compute the singular values exceeding
%  0.1 by applying deflation method.
%  [U,S,V] = svt(A,'lambda',0.1) 
%
%  Singular value thresholding, only compute the singular values exceeding
%  0.1 by applying succession method.
%  [U,S,V] = svt(A,'lambda',0.1,'method','succession')
%
%  Singular value decomposition for first 15 singular values, input is a
%  function handle, and the dimension of the original data matrix is
%  1000-by-1000.
%  [U,S,V] = svt(Afun,'k',15,'m',1000,'n',1000)
%
%  Singular value thresholding, only compute the singular values exceeding
%  0.1 by applying deflation method. Input is a function handle, and the
%  dimension of the original data matrix is 1000-by-1000.
%  [U,S,V] = svt(Afun,'lambda',0.1,'m',1000,'n',1000)
%
%  Singular value thresholding, only compute the singular values exceeding
%  0.1 by applying succession method. Input is a function handle, and the
%  dimension of the original data matrix is 1000-by-1000.
%  [U,S,V] = svt(Afun,'lambda',0.1,'m',1000,'n',1000,'method','succession')
%
%  More details can be seen in demo_svt.m
%
% COPYRIGHT: North Carolina State University 
% AUTHOR: Cai Li, Hua Zhou
% Email: cli9@ncsu.edu, hua_zhou@ncsu.edu

% Parse input
argin = inputParser;
argin.addRequired('A');
argin.addParamValue('lambda',NaN);
argin.addParamValue('k',6,@(x) x>0);
argin.addParamValue('incre',5,@(x) x>0);
argin.addParamValue('m',NaN);
argin.addParamValue('n',NaN);
argin.addParamValue('tol',1e-10,@(x) x>0);
argin.addParamValue('maxit',300,@(x) x>0);
argin.addParamValue('method','deflation',@(x) strcmp(x,'deflation') ...
                   || strcmp(x,'succession'))
argin.parse(A,varargin{:});

k = argin.Results.k;
incre = argin.Results.incre;
lambda = argin.Results.lambda;
tol = argin.Results.tol/sqrt(2);
maxit = argin.Results.maxit;
method = argin.Results.method;
if (strcmpi(method,'deflation'))
    def = 1;
else
    def = 0;
end

opts.tol = tol;
opts.maxit = maxit;

% Check input A type
if isnumeric(A) % If input A is a matrix
    if isnan(lambda) % Call svds directly for decomposition purpose 
        if nargout<=1
            opts.tol = opts.tol*sqrt(2);
            U = diag(svds(A,k,'L',opts));
        else
            opts.tol = opts.tol*sqrt(2);
            [U,S,V,flag] = svds(A,k,'L',opts);
        end
        return
    end
    [m,n] = size(A); % Obtain the dimension from the matrix
else % If input A is a function handle
    m = int32(argin.Results.m); % User need to input the dimension
    n = int32(argin.Results.n); % Dimension should be a positive integer
    ism = isa(m,'integer') && (m>0);
    isn = isa(n,'integer') && (n>0);
    if ~(ism && isn) % Check validity of input dimensions 
        error('Dimension must be a positive integer.');
    end
end

iter = min(m,n); % Set maximum succession number

% Check validity of k 
if (k>iter)
    k = iter;
    warning('K is out of dimension, reset to maximum value.');
end

% Check validity of lambda
if ~(isnan(lambda))
    if (isinf(lambda)) % Fast return for inf lambda
        U = []; S = []; V = []; flag = 0;
        return
    else
        if ~(lambda>=0)
            error('Lambda must be a nonnegative value.');
        end
    end
else
    if isnumeric(A)
        if nnz(A) == 0 % Fast return for zeros(m,n)
            if nargout<=1
                U = zeros(k,k);
            else
                U = eye(m,k);
                S = zeros(k,k);
                V = eye(n,k);
                flag = 0;
            end
            return
        end
    end
end

% Check validity of increment
if ~(isnan(lambda))
    if (incre>iter-k)
        incre = iter-k;
        warning('Incre is out of dimension; reset to maximum value.');
    end
end

w = [];  % Keep eigenvectors
e = [];  % Keep eigenvalues 
opts.issym = 1; % [zeros(n,n),A';zeros(m,m),A] is symmetric

% Main loop for computing singular values sequentially
while iter>0 
    [eigvecs,eigvals,eflag] = eigs(@matvec,double(m+n),double(k),'la',opts);
    eigvals = diag(eigvals); % Eigs output sorted
    if ~(isnan(lambda)) % For thresholding
        if (eflag) % Avoid non convergence situation ruin the deflation
            warning('eflag is %d, refresh with warm start.',eflag);
            k = max(length(e)-1,6); % Shift to avoid non convergence
            def = 0;
            [eigvecs,eigvals,eflag] = eigs(@matvec,double(m+n), ...
                double(k), 'la',opts);
            while(eflag)
                warning('eflag is %d, refresh with warm start.',eflag);
                k = max(k-1,1);
                [eigvecs,eigvals,eflag] = eigs(@matvec,double(m+n), ...
                double(k), 'la',opts);
            end
            eigvals = diag(eigvals); 
        end
    end
    if ~(isnan(lambda)) 
        i = find(eigvals<=lambda,1); % Thresholding
        if ~isempty(i) % Threshold found
            if def
                w = [w,eigvecs(:,1:i-1)]; %#ok<*AGROW>
                e = [e;eigvals(1:i-1)];
                if isempty(e) % Fast return for lambda>=max(eigvals)
                    U = []; S = []; V = []; flag = 0;
                    return
                end
            else
                w = eigvecs(:,1:i-1);
                e = eigvals(1:i-1);
                if isempty(e) % Fast return for lambda>=max(eigvals)
                    U = []; S = []; V = []; flag = 0;
                    return
                end
            end
            break
        end
    end
    if def % Deflation method adds the results every succession
        w = [w,eigvecs];
        e = [e;eigvals];
    else   % Succession method overwrites the results every succession 
        w = eigvecs;
        e = eigvals;
    end

    if (isnan(lambda)) % For decomposition
        break
    end
    
    iter = iter-k;
    if def % Increment based on deflation method
        k = min(incre,iter);
    else   % Increment based on succession method
        k = min(k+incre,iter);
    end
    incre = incre*2; % Non-constant increment
end

% Subfunction for function handle of eigs 
function mv = matvec(v)
    if isnumeric(A) % Matrix case
        mv = [(v(n+1:end)'*A)';A*v(1:n)]; % Avoid transpose of big matrix
    else  % Function handle case
        At = feval(A,v(n+1:end),true); 
        Af = feval(A,v(1:n),false);
        mv = [At;Af];
    end 
    
    if def % Deflation method
        if isempty(e)
            return;
        else
            mv = mv-w*(e.*(v'*w)'); % Deflation for eigen problem
        end
    end
end

% Output the results
if nargout<=1 % Output the singular values directly
    U = diag(e);
else
    S = diag(e); % Positive eigenvalues are singular values
    w = w*sqrt(2);
    
    if eflag % Tolerance to determine the small singular values
        dtol = max(m,n)*max(e)*sqrt(eps);
    else
        dtol = max(m,n)*max(e)*eps; 
    end

    j = find(abs(e)<=dtol,1); % Find the index of small singular values
    if isempty(j) 
        U = w(n+1:end,:); % Extract left singular vectors
        V = w(1:n,:); % Extract right singular vectors 
    else % Orthogonalize and extract the singular vectors
        U = [w(n+1:end,1:j-1),orth(w(n+1:end,j:end))]; 
        V = [w(1:n,1:j-1),orth(w(1:n,j:end))];
    end
    flag = eflag;
end

end