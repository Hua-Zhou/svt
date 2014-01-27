function[U,S,V,flag] = defsvt(A,varargin)
% Singular value thresholding based on deflation method
%
% USAGE:
%   [U,S,V,flag] = defsvt(A,'PARAM1',val1,'PARAM2',val2...)
%   S = defsvt(A,'PARAM1',val1,'PARAM2',val2...)
%
% INPUT:
%   A - m-by-n matrix or a function handle provided by user
%   lambda - threshold (default: NaN)
%   k - number of singular values to try (default: 6)
%   incre - increment of try after first k attemp (default: 3)
%   m - dimension of row
%   n - dimension of column
%   tol - eigs convergence tolerance (default: eps)
%   maxit - maximum number of eigs iterations (default: 300)
%
% OUTPUT:
%   U - left singular vectors
%   S - diagonal matrix of thresholded singular values
%   V - right singular vectors
%   flag - if 0, iterative eigs converged; 1, eigs not converged
%
% COPYRIGHT: North Carolina State University
% AUTHOR: Cai Li, Hua Zhou 
% Email: cli9@ncsu.edu

% Parse input
argin = inputParser;
argin.addRequired('A');
argin.addParamValue('lambda',NaN);
argin.addParamValue('k',6,@(x) x>0);
argin.addParamValue('incre',3,@(x) x>0);
argin.addParamValue('m',NaN);
argin.addParamValue('n',NaN);
argin.addParamValue('tol',eps,@(x) x>0);
argin.addParamValue('maxit',300,@(x) x>0);
argin.parse(A,varargin{:});

k = argin.Results.k;
incre = argin.Results.incre;
lambda = argin.Results.lambda;
tol = argin.Results.tol;
maxit = argin.Results.maxit;

opts.tol = tol;
opts.maxit = maxit;

% Check input A type
if isnumeric(A) % If input A is a matrix
    if isnan(lambda) % Call svds directly for non-thresholding matrix input 
        if nargout<=1
            U = diag(svds(A,k,'L',opts));
        else
            [U,S,V,flag] = svds(A,k,'L',opts);
        end
        return
    end
    [m,n] = size(A); % Obtain the dimension from the matrix
else % If input A is a function handle
    m = int8(argin.Results.m); % User need to input the dimension
    n = int8(argin.Results.n);
    ism = isa(m,'integer') && (m>0);
    isn = isa(n,'integer') && (n>0);
    if ~(ism && isn) % Check validity of input dimensions 
        error('Dimension must be a positive integer.');
    end
end

iter = min(m,n); % Set maximum iteration number

% Check validity of k 
if (k>iter)
    k = iter;
    warning('K is out of dimension, reseted to maximum value.');
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

% Check validity of incre
if ~(isnan(lambda))
    if (incre>iter-k)
        incre = iter-k;
        warning('Incre is out of dimension, reseted to maximum value.');
    end
end

w = [];  % Keep eigenvectors
e = [];  % Keep eigenvalues 
opts.issym = 1; % [zeros(n,n),A';zeros(m,m),A] is symmetric

% Main loop for computing singular values sequentially
while iter>0
    [eigvecs,eigvals,eflag] = eigs(@matvec,double(m+n),double(k),'la',opts); % Sort eigs output
    eigvals = diag(eigvals);
    if ~(isnan(lambda))
        i = find(eigvals<=lambda,1); % Thresholding
        if ~isempty(i)
            w = [w,eigvecs(:,1:i-1)];
            e = [e;eigvals(1:i-1)];
            break
        end
    end
    w = [w,eigvecs];
    e = [e;eigvals];

    if (isnan(lambda)) % Non-thresholding
        break
    end
    
    iter = iter-k;
    k = min(incre,iter);    
end

% Subfunction for function handle of eigs 
function mv = matvec(v)
    if isnumeric(A)
        mv = [A'*v(n+1:end);A*v(1:n)];
    else
        At = feval(A,v(n+1:end),true);
        Af = feval(A,v(1:n),false);
        mv = [At;Af];
    end 
    
    if isempty(e)
        return;
    else
        mv = mv-w*(e.*(w'*v)); % Deflation for eigen problem
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