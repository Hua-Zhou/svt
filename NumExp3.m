%% comparison with svds 
% singular value decomposition on structured matrix
clear;

% simulation settings 
rep = 11; 
seed = 2014;

% load testing matrices
load('bfwb398.mat');              % 398-by-398
load('rdb800l.mat');              % 800-by-800
load('tols1090.mat');             % 1090-by-1090
load('mhd4800b.mat');             % 4800-by-4800
load('cryg10000.mat');            % 10000-by-10000

% keeper of results 
records1 = zeros(rep-1,5);
records2 = zeros(rep-1,5);

% use Sim function to conduct comparison
records = Sim(bfwb398,'choice','stru-svds','rep',rep,'seed',seed);
records1(:,1) = records(:,1); % collect svt results
records2(:,1) = records(:,2); % collect objective function results

records = Sim(rdb800l,'choice','stru-svds','rep',rep,'seed',seed);
records1(:,2) = records(:,1);
records2(:,2) = records(:,2);

records = Sim(tols1090,'choice','stru-svds','rep',rep,'seed',seed);
records1(:,3) = records(:,1);
records2(:,3) = records(:,2);

records = Sim(mhd4800b,'choice','stru-svds','rep',rep,'seed',seed);
records1(:,4) = records(:,1);
records2(:,4) = records(:,2);

records = Sim(cryg10000,'choice','stru-svds','rep',rep,'seed',seed);
records1(:,5) = records(:,1);
records2(:,5) = records(:,2);

% save data
save('ne3_svt.mat','records1');
save('ne3_obj.mat','records2');

% %% boxplot
% D(:,1:2:2*size(records1,2)) = log(records1);
% D(:,2:2:2*size(records2,2)) = log(records2);
% 
% boxplot(D,'factorgap',10,'color','rk');
% set(gca,'xtick',1.9:3.9:50);
% set(gca,'xticklabel',...
%     {'bfwb398','rdb800l','tols1090','mhd4800b','cryg10000'});
% 
% ylabel('log(Run Time)','fontsize',20);
% legend(findobj(gca,'Tag','Box'),'svds','stru-svt','Location','northwest');
% set(gca,'fontsize',20);