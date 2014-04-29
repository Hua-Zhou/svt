%% application to matrix completion problems
% singular value thresholding on structured matrix
clear;

% simulation settings 
rep = 10; 
seed = 2014;
num = 10; % number of points in solution path

% keeper of results 
records1 = zeros(rep,5);
records2 = zeros(rep,5);
records3 = zeros(rep,5);

% use Sim function to conduct comparison
records = Sim_MatrixCompletion(500,500,'rep',rep,'seed',seed,'num',num);
records1(:,1) = records(:,1); % collect stru_svt results
records2(:,1) = records(:,2); % collect non_stru_svt results
records3(:,1) = records(:,3); % collect full svd results

records = Sim_MatrixCompletion(1000,1000,'rep',rep,'seed',seed,'num',num);
records1(:,2) = records(:,1);
records2(:,2) = records(:,2);
records3(:,2) = records(:,3);

records = Sim_MatrixCompletion(1500,1500,'rep',rep,'seed',seed,'num',num);
records1(:,3) = records(:,1);
records2(:,3) = records(:,2);
records3(:,3) = records(:,3);

records = Sim_MatrixCompletion(2000,2000,'rep',rep,'seed',seed,'num',num);
records1(:,4) = records(:,1);
records2(:,4) = records(:,2);
records3(:,4) = records(:,3);

records = Sim_MatrixCompletion(2500,2500,'rep',rep,'seed',seed,'num',num);
records1(:,5) = records(:,1);
records2(:,5) = records(:,2);
records3(:,5) = records(:,3);

% save data
save('ne9_svt.mat','records1');
save('ne9_non_stru.mat','records2');
save('ne_9_full_svd.mat','records3');

% %% boxplot
% D(:,1:3:3*size(records1,2)) = log(records1);
% D(:,2:3:3*size(records2,2)) = log(records2);
% D(:,3:3:3*size(records3,2)) = log(records3);
% 
% boxplot(D,'factorgap',10,'color','rkg');
% set(gca,'xtick',1.9:3.9:50);
% set(gca,'xticklabel',...
%     {'500','1000','1500','2000','2500'});
% 
% ylabel('log(Run Time)','fontsize',20);
% legend(findobj(gca,'Tag','Box'),'full svd','non-stru','stru-svt',...
% 'Location','northwest');
% set(gca,'fontsize',20);