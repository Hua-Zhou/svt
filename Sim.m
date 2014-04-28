function [records] = Sim(matrix_name,varargin)
% Parse input
argin = inputParser;
argin.addRequired('matrix_name');
argin.addParamValue('choice','svds',@(x) ischar(x));
argin.addParamValue('rep',5,@(x) x>0);
argin.addParamValue('seed',2014,@(x) x>0);
argin.addParamValue('rank',10,@(x) x>0);
argin.addParamValue('kth',50,@(x) x>0);
argin.addParamValue('load',false,@(x) islogical(x));
argin.parse(matrix_name,varargin{:});

choice = argin.Results.choice;
rep = argin.Results.rep;
seed = argin.Results.seed;
rank = argin.Results.rank;
kth = argin.Results.kth;
load = argin.Results.load;

s = RandStream('mt19937ar','Seed',seed);      % Reproducible
RandStream.setGlobalStream(s);                            

if (load)
    data = load('matrix_name'); % not work
    %data = Mat.Problem.A;           
else
    data = matrix_name;
end

% Basic feature about the sparse matrix
fprintf('The size of the matrix is %d by %d\n',size(data,1),size(data,2) );
fprintf('The sparsity of the matrix is %d\n', ...
    1-nnz(data)/(size(data,1)*size(data,2)));

records = zeros(2,rep);            % Store the run time

switch choice 
    case 'svd' %vs svd with sparse matrix
        fprintf('Compared with svd:\n');

        for j = 1:rep
           tic;
           [u,s,v] = svt(data);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(data-u*s*v','fro')/norm(data,'fro'));

           tic;
           [su,ss,sv] = svd(full(data));
           records(2,j) = toc;
           fprintf('Accuracy of svd: %d\n', ...
               norm(data-su(:,1:6)*ss(1:6,1:6)*sv(:,1:6)',...
               'fro')/norm(data,'fro'));
        end
        
    case 'svds' %vs svds with sparse matrix
        fprintf('Compared with svds:\n');

        for j = 1:rep
           tic;
           [u,s,v] = svt(data);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(data-u*s*v','fro')/norm(data,'fro'));

           tic;
           [su,ss,sv] = svds(data,6);
           records(2,j) = toc;
           fprintf('Accuracy of svds: %d\n', ...
               norm(data-su*ss*sv','fro')/norm(data,'fro'));
        end

    case 'stru-svds'    % vs svds with structrue matrix
        fprintf('Compared with svds for structured matrex:\n');
        m = size(data,1);
        n = size(data,2);
        l = randn(m,rank);  
        r = randn(n,rank);
        lr = l*r';           % Generation of low rank
        MA = data;
        new_data = MA + lr;    % Sparse + low rank
        
%         opts.tol = eps;
%         opts.maxit = 300;

        for j = 1:rep
           tic;
           [u,s,v] = svt(@MAtimesVec,'m',m,'n',n);
           records(1,j) = toc;
           fprintf('Accuracy of svt exploiting structure: %d\n', ...
               norm(new_data-u*s*v','fro')/norm(new_data,'fro'));

           tic;
           [su,ss,sv] = svds(new_data,6);
           records(2,j) = toc;
           fprintf('Accuracy of svds: %d\n', ...
               norm(new_data-su*ss*sv','fro')/norm(new_data,'fro'));
        end

    
    case 'svt-svd' % vs svd for svt with sparse matrix
        fprintf('Compared with svd for sparse matrix svt:\n');
        lambdas = svds(data,kth);
        tau = lambdas(kth);
        fprintf('lambda is: %d\n',tau);
        for j = 1:rep
           tic;
           [u,s,v] = svt(data,'lambda',tau);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(data-u*s*v','fro')/norm(data,'fro'));

           tic;
           [fu,fs,fv] = svd(full(data));
           fs = diag(fs);
           i = find(fs<=tau);
           fu = fu(:,1:i-1);
           fs = fs(1:i-1);
           fv = fv(:,1:i-1);
           fs = diag(fs);
           records(2,j) = toc;
           fprintf('Accuracy of full svd for svt: %d\n', ...
               norm(data-fu*fs*fv','fro')/norm(data,'fro'));
        end
    
    case 'svt-svd-stru' % vs svd for svt with structured matrix
        fprintf('Compared with svd for structured matrix svt:\n');
        m = size(data,1);
        n = size(data,2);
        l = randn(m,rank);  
        r = randn(n,rank);
        lr = l*r';           % Generation of low rank
        MA = data;
        new_data = MA + lr;    % Sparse + low rank
        
        lambdas = svds(new_data,kth);
        tau = lambdas(kth);
        fprintf('lambda is: %d\n',tau);
        
        for j = 1:rep
           tic;
           [u,s,v] = svt(@MAtimesVec,'lambda',tau,'m',m,'n',n);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(new_data-u*s*v','fro')/norm(new_data,'fro'));

           tic;
           [fu,fs,fv] = svd(full(new_data));
           fs = diag(fs);
           i = find(fs<=tau);
           fu = fu(:,1:i-1);
           fs = fs(1:i-1);
           fv = fv(:,1:i-1);
           fs = diag(fs);
           records(2,j) = toc;
           fprintf('Accuracy of full svd for svt: %d\n', ...
               norm(new_data-fu*fs*fv','fro')/norm(new_data,'fro'));
        end    
        
    case 'without-stru' % vs svt for not exploiting structure
        fprintf('Compared with svt for not exploiting structure:\n');
        
        m = size(data,1);
        n = size(data,2);
        l = randn(m,rank);  
        r = randn(n,rank);
        lr = l*r';           % Generation of low rank
        MA = data;
        new_data = MA + lr;    % Sparse + low rank
        
        lambdas = svds(new_data,kth);
        tau = lambdas(kth);
        fprintf('lambda is: %d\n',tau);

        for j = 1:rep
           tic;
           [u,s,v] = svt(@MAtimesVec,'lambda',tau,'m',m,'n',n);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(new_data-u*s*v','fro')/norm(new_data,'fro'));

           tic;
           [su,ss,sv] = svt(new_data,'lambda',tau);
           records(2,j) = toc;
           fprintf('Accuracy of svt not exploiting structure: %d\n', ...
               norm(new_data-su*ss*sv','fro')/norm(new_data,'fro'));
        end   
    
    case 'succession' %vs succession for svt with sparse matrix
        fprintf('Compared with succession:\n');
        lambdas = svds(data,kth);
        tau = lambdas(kth);
        fprintf('lambda is: %d\n',tau);
        for j = 1:rep
           tic;
           [u,s,v] = svt(data,'lambda',tau);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(data-u*s*v','fro')/norm(data,'fro'));

           tic;
           [su,ss,sv] = svt(data,'lambda',tau,'method','succession');
           records(2,j) = toc;
           fprintf('Accuracy of succession svt: %d\n', ...
               norm(data-su*ss*sv','fro')/norm(data,'fro'));
        end
     
    case 'stru-succession'    % vs succession with structrued matrix
        fprintf('Compared with stru-succession:\n');
        
        m = size(data,1);
        n = size(data,2);
        l = randn(m,rank);  
        r = randn(n,rank);
        lr = l*r';           % Generation of low rank
        MA = data;
        new_data = MA + lr;    % Sparse + low rank
        
        lambdas = svds(new_data,kth);
        tau = lambdas(kth);
        fprintf('lambda is: %d\n',tau);

        for j = 1:rep
           tic;
           [u,s,v] =svt(@MAtimesVec,'lambda',tau,'m',m,'n',n);
           records(1,j) = toc;
           fprintf('Accuracy of svt: %d\n', ...
               norm(new_data-u*s*v','fro')/norm(new_data,'fro'));

           tic;
           [su,ss,sv] = svt(@MAtimesVec,'lambda',tau,'m',m,'n',n, ...
               'method','succession');
           records(2,j) = toc;
           fprintf('Accuracy of succession svt: %d\n', ...
               norm(new_data-su*ss*sv','fro')/norm(new_data,'fro'));
        end

end
    
records = records(:,2:rep)';

fprintf('The mean of run time of svt is %d\n',mean(records(:,1)));
fprintf('The mean of run time of objective function is %d\n', ...
    mean(records(:,2)));
fprintf('The se of run time of svt is %d\n',std(records(:,1))/sqrt(rep-1));
fprintf('The se of run time of objective function is %d\n', ...
    std(records(:,2))/sqrt(rep-1));
fprintf('\n');


function MAvec = MAtimesVec(vec, trans) % Subfunction input by user 
    if trans
       MAvec = (vec'*MA)' + r*(vec'*l)';
    else
       MAvec = MA*vec + l*(r'*vec);
    end    
end


end