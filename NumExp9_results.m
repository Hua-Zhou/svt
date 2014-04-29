%% boxplot and summary
clear;

load('ne9_svt.mat');
load('ne9_non_stru.mat');
load('ne_9_full_svd.mat');

% mean
fprintf('The mean of run time of stru-svt is \n');
display(mean(records1));
fprintf('The mean of run time of non-stru-svt is \n');
display(mean(records2));
fprintf('The mean of run time of full svd is \n');
display(mean(records3));

% standard error
fprintf('The se of run time of svt is \n');
display(std(records1/sqrt(11-1)));
fprintf('The se of run time of non-stru-svt is \n');
display(std(records2/sqrt(11-1)));
fprintf('The se of run time of full svd is \n');
display(std(records3/sqrt(11-1)));

D(:,1:3:3*size(records1,2)) = log(records1);
D(:,2:3:3*size(records2,2)) = log(records2);
D(:,3:3:3*size(records3,2)) = log(records3);

boxplot(D,'factorgap',10,'color','rkg');
set(gca,'xtick',1.9:3.9:50); 
set(gca,'xticklabel',...
    {'500','1000','1500','2000','2500'});

ylabel('log(Run Time)','fontsize',20);
legend(findobj(gca,'Tag','Box'),'full svd','non-stru','stru-svt',...
'Location','northwest');
set(gca,'fontsize',20);