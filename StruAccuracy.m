function[] = StruAccuracy
% Reset random seed
s = RandStream('mt19937ar','Seed',2014);
RandStream.setGlobalStream(s);

clear;

data = load('bfwb398.mat');
mat = data.bfwb398;

% Generation of structured matrix (sparse plus low rank)
m = size(mat,1);
n = size(mat,2);
l = randn(m,20);  
r = randn(n,20);
lr = l*r';                
smat = mat + lr;

[u,s,v] = svt(@MAtimesVec,'lambda',1.849154e-05,'m',m,'n',n,'k',10);
fprintf('svt_de: %d\n',norm(smat-u*s*v','fro')/norm(smat,'fro'));
[su,ss,sv] = svds(smat,66);
fprintf('svds: %d\n',norm(smat-su*ss*sv','fro')/norm(smat,'fro'));
[eu,es,ev] = svt(@MAtimesVec,'lambda',1.849154e-05,'m',m,'n',n,'k',10,...
    'method','succession');
fprintf('svt: %d\n',norm(smat-eu*es*ev','fro')/norm(smat,'fro'));
fprintf('*************************\n');


clear;
data = load('rdb800l.mat');
mat = data.rdb800l;

% Generation of structured matrix (sparse plus low rank)
m = size(mat,1);
n = size(mat,2);
l = randn(m,20);  
r = randn(n,20);
lr = l*r';                
smat = mat + lr;

[u,s,v] = svt(@MAtimesVec,'lambda',2.758302e+01,'m',m,'n',n,'k',10);
fprintf('svt_de: %d\n',norm(smat-u*s*v','fro')/norm(smat,'fro'));
[su,ss,sv] = svds(smat,65);
fprintf('svds: %d\n',norm(smat-su*ss*sv','fro')/norm(smat,'fro'));
[eu,es,ev] = svt(@MAtimesVec,'lambda',2.758302e+01,'m',m,'n',n,'k',10,...
    'method','succession');
fprintf('svt: %d\n',norm(smat-eu*es*ev','fro')/norm(smat,'fro'));
fprintf('*************************\n');

clear;
data = load('tols1090.mat');
mat = data.tols1090;

% Generation of structured matrix (sparse plus low rank)
m = size(mat,1);
n = size(mat,2);
l = randn(m,20);  
r = randn(n,20);
lr = l*r';                
smat = mat + lr;

[u,s,v] = svt(@MAtimesVec,'lambda',1.115136e+06,'m',m,'n',n,'k',10);
fprintf('svt_de: %d\n',norm(smat-u*s*v','fro')/norm(smat,'fro'));
[su,ss,sv] = svds(smat,50);
fprintf('svds: %d\n',norm(smat-su*ss*sv','fro')/norm(smat,'fro'));
[eu,es,ev] = svt(@MAtimesVec,'lambda',1.115136e+06,'m',m,'n',n,'k',10,...
    'method','succession');
fprintf('svt: %d\n',norm(smat-eu*es*ev','fro')/norm(smat,'fro'));
fprintf('*************************\n');


% % clear;
% data = load('mhd4800b.mat');
% mat = data.mhd4800b;
% 
% % Generation of structured matrix (sparse plus low rank)
% m = size(mat,1);
% n = size(mat,2);
% l = randn(m,20);  
% r = randn(n,20);
% lr = l*r';                
% smat = mat + lr;
% 
% [u,s,v] = svt(@MAtimesVec,'lambda',9.302018e-02,'m',m,'n',n,'k',10);
% fprintf('svt_de: %d\n',norm(smat-u*s*v','fro')/norm(smat,'fro'));
% [su,ss,sv] = svds(smat,10);
% fprintf('svds: %d\n',norm(smat-su*ss*sv','fro')/norm(smat,'fro'));
% [eu,es,ev] = svt(@MAtimesVec,'lambda',9.302018e-02,'m',m,'n',n,'k',10,...
% 'method','succession');
% fprintf('svt: %d\n',norm(smat-eu*es*ev','fro')/norm(smat,'fro'));
% fprintf('*************************\n');

% data = load('cryg10000.mat');
% cryg10000 = data.cryg10000;
% [u,s,v] = svt(cryg10000,'lambda',2.199523e+04,'k',10);
% fprintf('svt_de: %d\n',norm(cryg10000-u*s*v','fro'));
% [su,ss,sv] = svds(cryg10000,50);
% fprintf('svds: %d\n',norm(cryg10000-su*ss*sv','fro'));
% [eu,es,ev] = svt(cryg10000,'lambda',2.199523e+04,'k',10,'method','deflation');
% fprintf('svt: %d\n',norm(cryg10000-eu*es*ev','fro'));
% fprintf('*************************\n');
% 
% data = load('af23560.mat');
% af23560 = data.af23560;
% [u,s,v] = svt(af23560,'lambda',3.636846e+02,'k',10);
% fprintf('svt: %d\n',norm(af23560-u*s*v','fro'));
% [su,ss,sv] = svds(af23560,50);
% fprintf('svds: %d\n',norm(af23560-su*ss*sv','fro'));
% [eu,es,ev] = svt(af23560,'lambda',3.636846e+02,'k',10,'method','deflation');
% fprintf('svt: %d\n',norm(af23560-eu*es*ev','fro'));
% fprintf('*************************\n');
 
% cond = load('cond-mat-2003');
% [u,s,v] = svt(condmat2003,'lambda',2.125883e+01,'k',50);
% fprintf('svt: %d\n',norm(condmat2003-u*s*v','fro'));
% [su,ss,sv] = svds(condmat2003,50);
% fprintf('svds: %d\n',norm(condmat2003-su*ss*sv','fro'));
% 
% c = load('c-66');
% [u,s,v] = svt(c66,'lambda',1.004267e+03,'k',50);
% fprintf('svt: %d\n',norm(c66-u*s*v','fro'));
% [su,ss,sv] = svds(c66,50);
% fprintf('svds: %d\n',norm(c66-su*ss*sv','fro'));

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