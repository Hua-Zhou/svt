function[] = GeneralComparison

s = RandStream('mt19937ar','Seed',2014);  % Reset random seed
RandStream.setGlobalStream(s);

p1 = 1000; % Dimension
p2 = 1000;
r = 30;     % Rank
s = 0.99;  % Sparsity
gridpts = 5;  % Grid points

l = randn(p1,r);
r = randn(r,p2);
lr = l * r;                   % Low rank matrix generation
display(rank(lr));

m = randsparse(p1,p2,s);      % Sparse matrix generation
% data = load('tols1090.mat'); 
% m = data.tols1090;
display(1-nnz(m)/numel(m));

mat = m + lr;                 % Data matrix generation

maxlambda = svds(mat,1);     
disp(maxlambda);

% lambdas = exp(log(maxlambda)/gridpts*(gridpts:-1:1));
lambdas = maxlambda-10*(gridpts:-1:1);

display('svt');
tic;
%profile on;
for i=1:gridpts
    [sU,ss,sV] = svt(mat,'lambda',lambdas(i),'method','succession');
    display(['Grid point ' num2str(i) ]);
end
%profile viewer;
toc;

display('structure_svt');
tic;
%profile on;
for i=1:gridpts
    [stU,sts,stV] = svt(@MAtimesVec,'m',p1,'n',p2,'lambda',lambdas(i),'method','succession');
    display(['Grid point ' num2str(i) ]);
end
%profile viewer;
toc;

display('full');
tic;
for i=1:gridpts
    [fU,fs,fV] = fsvt(mat,lambdas(i)); 
    display(['Grid point ' num2str(i) ]);
end
toc;

% display('full');
% tic;
% [fU2,fs2,fV2] = fsvt(full(m),lambdas(4));
% toc;

% tic;
% [fU2,fs2,fV2] = svd(full(m));
% toc;

% display('structure_svt');
% tic;
% profile on;
% [stU2,sts2,stV2] = svt(@MAtimesVec,'m',p1,'n',p2,'lambda',lambdas(4),'method','succession');
% profile viewer;
% toc;
% 
% display('svt');
% tic;
% %profile on;
% [sU20,ss20,sV20] = svt(m,'lambda',lambdas(4));
% %profile viewer;
% toc;
% 
% display('svt');
% tic;
% %profile on;
% [sU2,ss2,sV2] = svt(m,'lambda',lambdas(4),'method','succession');
% %profile viewer;
% toc;

function MAvec = MAtimesVec(vec, trans)
   if trans
     MAvec = (vec'*m)' + r'*(vec'*l)';
   else
     MAvec = m*vec + l*(r*vec);
   end
end

end
