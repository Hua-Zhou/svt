%% boxplot and summary
clear;

load('ne5_svt.mat');
load('ne5_obj.mat');

% mean
fprintf('The mean of run time of svt is \n');
display(mean(records1));
fprintf('The mean of run time of obj is \n');
display(mean(records2));

% standard error
fprintf('The se of run time of svt is \n');
display(std(records1/sqrt(11-1)));
fprintf('The se of run time of obj is \n');
display(std(records2/sqrt(11-1)));

D(:,1:2:2*size(records1,2)) = log(records1);
D(:,2:2:2*size(records2,2)) = log(records2);

boxplot(D,'factorgap',10,'color','rk');
set(gca,'xtick',1.9:3.9:50);
set(gca,'xticklabel',...
    {'bfwb398','rdb800l','tols1090','mhd4800b','cryg10000'});

ylabel('log(Run Time)','fontsize',20);
legend(findobj(gca,'Tag','Box'),'svd','stru-svt','Location','northwest');
set(gca,'fontsize',20);