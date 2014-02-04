% Testing, lambda is the 50th singular value
load('test.mat');

% size == 23560
profile on;
[u1,s1,v1]=defsvt(af23560,'lambda',3.636846e+02,'k',10,'incre',5);
profile off;
profile viewer;

% size == 31163
% profile on;
% [u2,s2,v2]=defsvt(condmat2003,'lambda',2.125883e+01,'k',10,'incre',5);
% profile off;
% profile viewer;