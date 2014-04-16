clear;

rep = 5; 
load('bfwb398.mat');
load('rdb800l.mat');
load('tols1090.mat');
load('mhd4800b.mat');
% records1 = zeros(rep-1,1);
% records2 = zeros(rep-1,1);

records_11 = Sim(bfwb398,'choice','svds','rep',rep,'seed',2014);
records_12 = Sim(bfwb398,'choice','stru-svds','rep',rep,'seed',2014);
records_13 = Sim(bfwb398,'choice','svt-svd','rep',rep,'seed',2014);
records_14 = Sim(bfwb398,'choice','svt-svd-stru','rep',rep,'seed',2014);
records_15 = Sim(bfwb398,'choice','without-stru','rep',rep,'seed',2014);
records_16 = Sim(bfwb398,'choice','succession','rep',rep,'seed',2014);
records_17 = Sim(bfwb398,'choice','stru-succession','rep',rep,'seed',2014);

records_21 = Sim(rdb800l,'choice','svds','rep',rep,'seed',2014);
records_22 = Sim(rdb800l,'choice','stru-svds','rep',rep,'seed',2014);
records_23 = Sim(rdb800l,'choice','svt-svd','rep',rep,'seed',2014);
records_24 = Sim(rdb800l,'choice','svt-svd-stru','rep',rep,'seed',2014);
records_25 = Sim(rdb800l,'choice','without-stru','rep',rep,'seed',2014);
records_26 = Sim(rdb800l,'choice','succession','rep',rep,'seed',2014);
records_27 = Sim(rdb800l,'choice','stru-succession','rep',rep,'seed',2014);

records_31 = Sim(tols1090,'choice','svds','rep',rep,'seed',2014);
records_32 = Sim(tols1090,'choice','stru-svds','rep',rep,'seed',2014);
records_33 = Sim(tols1090,'choice','svt-svd','rep',rep,'seed',2014);
records_34 = Sim(tols1090,'choice','svt-svd-stru','rep',rep,'seed',2014);
records_35 = Sim(tols1090,'choice','without-stru','rep',rep,'seed',2014);
records_36 = Sim(tols1090,'choice','succession','rep',rep,'seed',2014);
records_37 = Sim(tols1090,'choice','stru-succession','rep',rep,'seed',2014);

records_41 = Sim(mhd4800b,'choice','svds','rep',rep,'seed',2014);
records_42 = Sim(mhd4800b,'choice','stru-svds','rep',rep,'seed',2014);
records_43 = Sim(mhd4800b,'choice','svt-svd','rep',rep,'seed',2014);
records_44 = Sim(mhd4800b,'choice','svt-svd-stru','rep',rep,'seed',2014);
records_45 = Sim(mhd4800b,'choice','without-stru','rep',rep,'seed',2014);
records_46 = Sim(mhd4800b,'choice','succession','rep',rep,'seed',2014);
records_47 = Sim(mhd4800b,'choice','stru-succession','rep',rep,'seed',2014);