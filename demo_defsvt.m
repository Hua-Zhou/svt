%% Singular value decomposition for sparse matrix

clear;
% reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% read in sparse matrices downloaded from The University of Florida Sparse Matrix Collection 
load('bfwb398.mat');

% load('cryg10000.mat');
% load('af23569.mat');

mat = bfwb398;
fprintf('Size of bfwb398: ');
disp(size(mat));
fprintf('Sparsity of bfwb398: ');
disp(1-nnz(mat)/(size(mat,1)*size(mat,2))) % sparsity

%%
% comparison with svds
fprintf('Run time of defsvt: ');
tic;
[u,s,v] = defsvt(mat,'k',10);
disp(toc);

fprintf('Run time of svds: ');
tic;
[su,ss,sv] = svds(mat,10);
disp(toc);

fprintf('Accuracy of solutions provied by defsvt');
disp(norm(mat-u*s*v','fro')/norm(mat,'fro')); % accuracy
fprintf('Accuracy of solutions provied by svds');
disp(norm(mat-su*ss*sv','fro')/norm(mat,'fro'));

%% Singular value decomposition for structured but non-sparse matrix (sparse plus low rank)

clear;
% reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% read in sparse matrices downloaded from The University of Florida Sparse Matrix Collection 
load('tols1090.mat');

mat = tols1090;
fprintf('Size of tols1090: ');
disp(size(mat));
fprintf('Sparsity of tols1090: ');
disp(1-nnz(mat)/(size(mat,1)*size(mat,2)))

%%
% generation of non-sparse matrix (sparse plus low rank)
m = size(mat,1);
n = size(mat,2);
p = randn(m,20);  
l = randn(n,20);
pl = p*l';                % generation of low rank matrix
smat = mat + pl;          % sparse + low rank
fprintf('Rank: ');
disp(rank(pl));

%%
% comparison with svds
fprintf('Run time of defsvt: ');
%tic;
%[u,s,v] = defsvt(@MAtimesVec,'m',m,'n',n,'k',10);
%disp(toc);

fprintf('Run time of svds: ');
tic;
[su,ss,sv] = svds(smat,10);
disp(toc);

%fprintf('Accuracy of solutions provied by defsvt');
%disp(norm(smat-u*s*v','fro')/norm(smat,'fro'));
fprintf('Accuracy of solutions provied by svds');
disp(norm(smat-su*ss*sv','fro')/norm(smat,'fro'));

%% Singular value thresholding for sparse matrix

clear;
% reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

%%
% read in sparse matrices downloaded from The University of Florida Sparse Matrix Collection 
load('mhd4800b.mat');
%load('rdb800l.mat');

%mat = rdb800l;
mat = mhd4800b;
fprintf('Size of rdb800l: ');
disp(size(mat));
fprintf('Sparsity of rdb800l: ');
disp(1-nnz(mat)/(size(mat,1)*size(mat,2)))

%%
% comparison with svd & iterative svds
fprintf('Run time of defsvt: ');
tic;
[u,s,v] = defsvt(mat,'lambda',4.152050e-02);
disp(toc);

fprintf('Run time of svd: ');
%tic;
%[su,ss,sv] = svd(full(mat));
%disp(toc);

fprintf('Run time of itersvt: ');
tic;
[iu,is,iv] = itersvt(mat,'lambda',4.152050e-02);
disp(toc);

% fprintf('Accuracy of solutions provied by defsvt');
% disp(norm(mat-u*s*v','fro')/norm(mat,'fro'));
% fprintf('Accuracy of solutions provied by svd');
% disp(norm(mat-su*ss*sv','fro')/norm(mat,'fro'));
% fprintf('Accuracy of solutions provied by itersvt');
% disp(norm(mat-iu*is*iv','fro')/norm(mat,'fro'));

%% Singular value thresholding for structured but non-sparse matrix (sparse plus low rank)

% clear;
% reset random seed
% s = RandStream('mt19937ar','Seed',2014);
% RandStream.setGlobalStream(s);
% 
% %
% read in sparse matrices downloaded from The University of Florida Sparse Matrix Collection 
% load('mhd4800b.mat');
% 
% mat = mhd4800b;
% fprintf('Size of mhd4800b: ');
% disp(size(mat));
% fprintf('Sparsity of mhd4800b: ');
% disp(1-nnz(mat)/(size(mat,1)*size(mat,2)))
% 
% %
% generation of non-sparse matrix (sparse plus low rank)
% m = size(mat,1);
% n = size(mat,2);
% p = randn(m,20);  
% l = randn(n,20);
% pl = p*l';                % generation of low rank matrix
% smat = mat + pl;          % sparse + low rank
% fprintf('Rank: ');
% disp(rank(pl));

%%
% comparison with iterative svds
% fprintf('Run time of defsvt: ');
% tic;
% [u,s,v] = defsvt(@MAtimesVec,'m',m,'n',n,'lambda',4.152050e-02);
% disp(toc);

% fprintf('Run time of itersvt: ');
% tic;
% [iu,is,iv] = itersvt(@MAtimesVec,'m',m,'n',n,'lambda',4.152050e-02);
% disp(toc);
% 
% fprintf('Accuracy of solutions provied by defsvt');
% disp(norm(smat-u*s*v','fro')/norm(smat,'fro'));
% fprintf('Accuracy of solutions provied by itersvt');
% disp(norm(smat-iu*is*iv','fro')/norm(smat,'fro'));





