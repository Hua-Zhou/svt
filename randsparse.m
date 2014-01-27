function mat = randsparse(m,n,s)
% m: number of rows
% n: number of columns
% s: expected sparsity

num = round(m*n*(1-s));
rows = randi(m,num,1);
cols = randi(n,num,1);
vec = randi(10,num,1);
mat = sparse(rows,cols,vec);

end