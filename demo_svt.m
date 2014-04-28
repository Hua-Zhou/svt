%%
% We need a functional environment so that the MAtimesVec subfunction can
% access variables in workspace
function[] = demo_svt

%% Singular value decomposition for sparse matrix

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
% Size of matrix
disp(size(mat));

%%
% Sparsity of matrix
disp(nnz(mat)/numel(mat)); 

%%
% Top 25 singular values/vectors by svt
tic;
[u,s,v] = svt(mat,'k',25);
toc;

%%
% Top 25 singular values/vectors by Matlab svds
tic;
[su,ss,sv] = svds(mat,25);
toc;

%%
% Full svd
tic;
[fu,fs,fv] = svd(full(mat));
toc;

%%
% Accuracy of solutions provided by svt
disp(norm(fu(:,1:25)*fs(1:25,1:25)*fv(:,1:25)'- u*s*v','fro')); 

%%
% Accuracy of solutions provided by svds
disp(norm(fu(:,1:25)*fs(1:25,1:25)*fv(:,1:25)'- su*ss*sv','fro')); 

%% Singular value decomposition for structured (sparse + low rank) matrix

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
L = randn(m,20);  
R = randn(n,20);
LR = L*R';                % generation of low rank matrix
smat = mat + LR;          % sparse + low rank

%%
% Top 25 singular values/vectors by svt. Function MAtimesVec is defined at
% end of this file.
tic;
[u,s,v] = svt(@MAtimesVec,'m',m,'n',n,'k',25);
toc;

%%
% Top 25 singular values/vectors by Matlab's svds
tic;
[su,ss,sv] = svds(smat,25);
toc;

%%
% Full svd
tic;
[fu,fs,fv] = svd(full(smat));
toc;

%%
% Accuracy of solutions provided by svt
disp(norm(fu(:,1:25)*fs(1:25,1:25)*fv(:,1:25)'- u*s*v','fro')); 

%%
% Accuracy of solutions provided by svds
disp(norm(fu(:,1:25)*fs(1:25,1:25)*fv(:,1:25)'- su*ss*sv','fro')); 

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
disp(nnz(mat)/numel(mat));

%%
% Find all singular values >= 0.1 by svt (deflation method)
tic;
[u,s,v] = svt(mat,'lambda',0.1);
toc;
display(size(s));

%%
% It's faster if we have a good guess of how many singular values above
% threshold
tic;
[~,ks,~] = svt(mat,'lambda',0.1,'k',45);
toc;
display(size(ks));

%%
% Find all singular values >= 0.1 by svt (succession method)
tic;
[iu,is,iv] = svt(mat,'lambda',0.1,'method','succession');
toc;
display(size(is));

%%
% Find all singular values >= 0.1 by full svd
fmat = full(mat);
tic;
[su,ss,sv] = svd(fmat);
dss = diag(ss);
i = find(dss<=0.1);
su = su(:,1:i-1);
dss = dss(1:i-1);
sv = sv(:,1:i-1);
ss = diag(dss);
toc;
display(size(ss));

%%
% Accuracy of solutions provided by svt deflation method
disp(norm(u*s*v'- su*ss*sv','fro')); 

%%
% Accuracy of solutions provided by svt succession method
disp(norm(iu*is*iv'- su*ss*sv','fro')); 

%% Singular value thresholding for structured (sparse + low rank) matrix 

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
L = randn(m,20);  
R = randn(n,20);
LR = L*R';                % generation of low rank matrix
smat = mat + LR;          % sparse + low rank

%%
% Find all singular values >= 0.2 by svt (deflation method). Function
% MAtimesVec is defined at end of this file.
tic;
[u,s,v] = svt(@MAtimesVec,'m',m,'n',n,'lambda',0.2);
toc;
display(size(s));

%%
% It's faster if we have a good guess of how many singular values above
% threshold
tic;
[~,ks,~] = svt(@MAtimesVec,'m',m,'n',n,'lambda',0.2,'k',45);
toc;
display(size(ks));

%%
% Find all singular values >= 0.2 by svt (succession method). Function
% MAtimesVec is defined at end of this file.
tic;
[iu,is,iv] = svt(@MAtimesVec,'m',m,'n',n,'lambda',0.2,...
'method','succession');
toc;
display(size(is));

%%
% Find all singular values >= 0.2 by full svd
fmat = full(smat);
tic;
[su,ss,sv] = svd(fmat);
dss = diag(ss);
i = find(dss<=0.2);
su = su(:,1:i-1);
dss = dss(1:i-1);
sv = sv(:,1:i-1);
ss = diag(dss);
toc;
display(size(ss));

%%
% Accuracy of solutions provided by svt deflation method
disp(norm(u*s*v'- su*ss*sv','fro')); 

%%
% Accuracy of solutions provided by svt succession method
disp(norm(iu*is*iv'- su*ss*sv','fro')); 

%%
% Subfunction for exploiting matrix structure of sparse plus low rank
function MAvec = MAtimesVec(vec, trans)

    if trans
       MAvec = (vec'*mat)' + R*(vec'*L)';
    else
       MAvec = mat*vec + L*(R'*vec);
    end
    
end

end




