function [u,s,v,f,data]=funhandle
%MA = [0,0,0,1,2,0;0,3,0,1,9,0;0,1,0,0,0,1;3,1,2,1,0,0];
%MA = [0,0,1,3,0,0,1;0,2,4,1,9,0,1;0,0,0,0,1,9,2;9,0,1,0,0,3,0;0,0,1,0,0,3,7];
%   west = load('west2021');
%   MA = west.Problem.A;
%   MA = MA(1:25,1:15);
u = randn(63,3);
l = randn(63,3);

MA = randsparse(63,63,0.9);
data = MA + u*l';

%[su,ss,sv] = svds(A);
[u,s,v,f] = defsvt(@MAtimesVec,'m',size(data,1),'n',size(data,2));

function MAvec = MAtimesVec(vec, varargin)
    argin = inputParser;
    argin.addRequired('vec');
    argin.addOptional('trans', false, @islogical);
    argin.parse(vec,varargin{:});

    trans = argin.Results.trans;

    if trans
       %MAvec = (vec'*MA)';
       MAvec = (vec'*MA)' + l*(vec'*u)';
    else
       %MAvec = MA*vec; 
       MAvec = MA*vec + u*(l'*vec);
    end
    
end





end