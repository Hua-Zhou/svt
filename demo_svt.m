function[] = demo_svt
%% Singular value decomposition for sparse matrix

clear;
% Reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% Read in sparse matrices downloaded from The University of Florida Sparse
% Matrix Collection
data = load('tols1090.mat');
mat = data.tols1090;

%%
% Size of matrix
disp(size(mat));

%%
% Sparsity of matrix
disp(nnz(mat)/numel(mat)); 

%%
% Run time of svt for top 10 singular value decomposition: 
tic;
[u,s,v] = svt(mat,'k',10);
disp(toc);

%%
% Run time of svds for top 10 singular value decomposition: 
tic;
[su,ss,sv] = svds(mat,10);
disp(toc);

%%
% Accuracy of solutions provided by svt
disp(norm(mat-u*s*v','fro')/norm(mat,'fro')); 

%%
% Accuracy of solutions provided by svds
disp(norm(mat-su*ss*sv','fro')/norm(mat,'fro'));

%% Singular value decomposition for structured matrix

clear;
% Reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% Read in sparse matrices downloaded from The University of Florida Sparse
% Matrix Collection
data = load('tols1090.mat');
mat = data.tols1090;

%%
% Generation of structured matrix (sparse plus low rank)
m = size(mat,1);
n = size(mat,2);
l = randn(m,20);  
r = randn(n,20);
lr = l*r';                % generation of low rank matrix
smat = mat + lr;          % sparse + low rank
disp(rank(lr));           % rank 

%%
% Run time of svt for top 10 singular value decomposition
tic;
[u,s,v] = svt(@MAtimesVec,'m',m,'n',n,'k',10);
disp(toc);

%%
% Run time of svds for top 10 singular value decomposition
tic;
[su,ss,sv] = svds(smat,10);
disp(toc);

%%
% Accuracy of solutions provided by svt
disp(norm(smat-u*s*v','fro')/norm(smat,'fro'));

%%
% Accuracy of solutions provided by svds
disp(norm(smat-su*ss*sv','fro')/norm(smat,'fro'));

%% Singular value thresholding for sparse matrix

clear;
% Reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% Read in sparse matrices downloaded from The University of Florida Sparse
% Matrix Collection
data = load('mhd4800b.mat');
mat = data.mhd4800b;

%%
% Size of matrix:
disp(size(mat));

%%
% Sparsity of matrix
disp(1-nnz(mat)/numel(mat));

%%
% Run time of svt for singular value thresholding by applying deflation
% method
tic;
[u,s,v] = svt(mat,'lambda',4.152050e-02);
disp(toc);

%%
% Run time of svd for full singular value thresholding 
fmat = full(mat);
tic;
[su,ss,sv] = svd(fmat);
dss = diag(ss);
i = find(dss<=4.152050e-02);
su = su(:,1:i-1);
dss = dss(1:i-1);
sv = sv(:,1:i-1);
ss = diag(dss);
disp(toc);

%%
% Run time of svt for singular value thresholding by applying succession
% method
tic;
[iu,is,iv] = svt(mat,'lambda',4.152050e-02,'method','succession');
disp(toc);

%%
% Accuracy of solutions provided by svt based on deflation method 
disp(norm(mat-u*s*v','fro')/norm(mat,'fro'));

%%
% Accuracy of solutions provided by svd  
disp(norm(fmat-su*ss*sv','fro')/norm(fmat,'fro'));

%%
% Accuracy of solutions provided by svt based on succession method
disp(norm(mat-iu*is*iv','fro')/norm(mat,'fro'));

%% Singular value thresholding for structured matrix 

clear;
% Reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% Read in sparse matrices downloaded from The University of Florida Sparse
% Matrix Collection
data = load('mhd4800b.mat'); 
mat = data.mhd4800b;

%%
% Generation of structured matrix (sparse plus low rank)
m = size(mat,1);
n = size(mat,2);
l = randn(m,20);  
r = randn(n,20);
lr = l*r';                % generation of low rank matrix
smat = mat + lr;          % sparse + low rank
disp(rank(lr));           % rank

%%
% Run time of svt for singular value thresholding by applying deflation
% method
tic;
[u,s,v] = svt(@MAtimesVec,'m',m,'n',n,'lambda',4.152050e-02);
disp(toc);

%%
% Run time of svd for full singular value thresholding 
fmat = full(smat);
tic;
[su,ss,sv] = svd(fmat);
dss = diag(ss);
i = find(dss<=4.152050e-02);
su = su(:,1:i-1);
dss = dss(1:i-1);
sv = sv(:,1:i-1);
ss = diag(dss);
disp(toc);

%%
% Run time of svt for singular value thresholding by applying deflation
% method
tic;
[iu,is,iv] = svt(@MAtimesVec,'m',m,'n',n,'lambda',4.152050e-02,...
'method','succession');
disp(toc);

%%
% Accuracy of solutions provided by svt based on deflation method 
disp(norm(smat-u*s*v','fro')/norm(smat,'fro'));

%%
% Accuracy of solutions provided by svd  
disp(norm(fmat-su*ss*sv','fro')/norm(fmat,'fro'));

%%
% Accuracy of solutions provided by svt based on succession method 
disp(norm(smat-iu*is*iv','fro')/norm(smat,'fro'));

%%
% Subfunction for utilizing matrix structure of sparse plus low rank
function MAvec = MAtimesVec(vec, varargin)
    argin = inputParser;
    argin.addRequired('vec');
    argin.addOptional('trans', false, @islogical);
    argin.parse(vec,varargin{:});

    trans = argin.Results.trans;

    if trans
       MAvec = (vec'*mat)' + r*(vec'*l)';
    else
       MAvec = mat*vec + l*(r'*vec);
    end
    
end

end




