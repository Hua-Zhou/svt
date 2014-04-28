function [records] = Sim_MatrixCompletion(p1,p2,varargin)
% Parse input
argin = inputParser;
argin.addRequired('p1');
argin.addRequired('p2');
argin.addParamValue('rep',5,@(x) x>0);
argin.addParamValue('seed',2014,@(x) x>0);
argin.addParamValue('rank',50,@(x) x>0);
argin.addParamValue('prop',0.95,@(x) x>0 && x<1);
argin.addParamValue('gridpts',50,@(x) x>0);
argin.addParamValue('num',5,@(x) x>0);
argin.parse(p1,p2,varargin{:});

rep = argin.Results.rep;
seed = argin.Results.seed;
r = argin.Results.rank;
missingprop = argin.Results.prop;
gridpts = argin.Results.gridpts;
num = argin.Results.num;

s = RandStream('mt19937ar','Seed',seed);      % Reproducible
RandStream.setGlobalStream(s);  

M = randn(p1,r) * randn(r,p2) + randn(p1,p2); % Plus random noise
display(rank(M));

missingidx = rand(p1,p2)<missingprop;

Mobs = M;
Mobs(missingidx) = nan;

records = zeros(3,rep);     % Keeper
% find max lambda for nuclear norm regularization
[~,stats] = MatrixCompletion_MM(Mobs,inf);
maxlambda = stats.maxlambda;
disp(maxlambda);

% Solution path generated
lambdas = exp(log(maxlambda)/gridpts*(gridpts:-1:1));

%% svt exploiting matrix structure
display('structure_svt');
% solution path by warm start
for j = 1:rep
    fprintf('**********Replication %d**********\n',j)
    tic;
    %profile on;
    for i=1:num
        if (i==1)
            Y0 = [];
        else
            Y0 = Z;
        end
        [Z,stats] = MatrixCompletion_MM(Mobs,lambdas(i),'Y0',Y0,'Display','off');

        if i >= 3
            display(['Grid point ' num2str(i) ', rank=' num2str(stats.rank)]);
        end
    end
    %profile viewer;
    records(1,j) = toc;
end

%% svt w/o exploiting matrix structure
display('svt');
tic;
%profile on;
for j = 1:rep
    fprintf('**********Replication %d**********\n',j)
    tic;
    %profile on;
    for i=1:num
        if (i==1)
            Y0 = [];
        else
            Y0 = Z;
        end
        [Z,stats] = MatrixCompletion_MM(Mobs,lambdas(i),'Y0',Y0,...
            'Display','off','method','svt');

        if i >= 3
            display(['Grid point ' num2str(i) ', rank=' num2str(stats.rank)]);
        end
    end
    %profile viewer;
    records(2,j) = toc;
end

%% full svt via svd
display('full_svt');
tic;
%profile on;
for j = 1:rep
    fprintf('**********Replication %d**********\n',j)
    tic;
    %profile on;
    for i=1:num
        if (i==1)
            Y0 = [];
        else
            Y0 = Z;
        end
        [Z,stats] = MatrixCompletion_MM(Mobs,lambdas(i),'Y0',Y0,...
            'Display','off','method','full');

        if i >= 3
            display(['Grid point ' num2str(i) ', rank=' num2str(stats.rank)]);
        end
    end
    %profile viewer;
    records(3,j) = toc;
end

records = records(:,1:rep)';  % Discard the 1st rep?
fprintf('**********Summary**********\n');
fprintf('The mean of run time of stru_svt is %d\n',mean(records(:,1)));
fprintf('The mean of run time of non_stru_svt is %d\n', mean(records(:,2)));
fprintf('The mean of run time of full svt is %d\n', mean(records(:,3)));

fprintf('The se of run time of stru_svt is %d\n',std(records(:,1))/sqrt(rep));
fprintf('The se of run time of non_stru_svt is %d\n', std(records(:,2))/sqrt(rep));
fprintf('The se of run time of full svt is %d\n', std(records(:,3))/sqrt(rep));
fprintf('\n');

end