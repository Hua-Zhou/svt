%% generate a random matrix with missing entries

clear;

s = RandStream('mt19937ar','Seed',2014);  % Reset random seed
RandStream.setGlobalStream(s);

p1 = 1500;
p2 = 1500;
r = 5;

M = randn(p1,r) * randn(r,p2);
display(rank(M));

missingprop = 0.75;
missingidx = rand(p1,p2)<missingprop;

Mobs = M;
Mobs(missingidx) = nan;

%% impute by SVT (Candes et al paper)
% 
% Z = matrix_impute(Mobs,'NoiseTol',.05,'Tau',500);
% Z = round(Z);
% missrate = nnz(Z(missingidx)~=M(missingidx))/nnz(missingidx);
% disp(missrate);

%% find max lambda for nuclear norm regularization

%[~,stats] = matrix_impute_Nesterov(Mobs,inf);
[~,stats] = matrix_impute_MM(Mobs,inf);
maxlambda = stats.maxlambda;
disp(maxlambda);

%%
gridpts = 11;
lambdas = exp(log(maxlambda)/gridpts*(gridpts:-1:1));

% solution path by warm start
tic;
for i=1:gridpts
    if (i==1)
        Y0 = [];
    else
        Y0 = Z;
    end
    %[Z,stats] = matrix_impute_Nesterov(Mobs,lambdas(i),'Y0',Y0,'Display','off');
    [Z,stats] = matrix_impute_MM(Mobs,lambdas(i),'Y0',Y0,'Display','off');
    
    if i >=3
        display(['Grid point ' num2str(i) ', rank=' num2str(stats.rank)]);
    end
end
toc;


