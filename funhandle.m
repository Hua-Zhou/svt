function [u,s,v,f,data] = funhandle

%MA = randsparse(3000,3000,0.9);

west = load('west2021');
MA = west.Problem.A;
u = randn(size(MA,1),100);
l = randn(size(MA,2),100);


data = MA + u*l';

[u,s,v,f] = defsvt(@MAtimesVec,'m',size(data,1),'n',size(data,2));

function MAvec = MAtimesVec(vec, varargin)
    argin = inputParser;
    argin.addRequired('vec');
    argin.addOptional('trans', false, @islogical);
    argin.parse(vec,varargin{:});

    trans = argin.Results.trans;

    if trans
       MAvec = (vec'*MA)' + l*(vec'*u)';
    else
       MAvec = MA*vec + u*(l'*vec);
    end
    
end





end