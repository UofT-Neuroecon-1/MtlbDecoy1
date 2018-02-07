function [m,v] = pascalstat(n,p)
%PASCALSTAT Mean and variance of the Pascal Distribution
%   [M,V] = PASCALSTAT(N,P) returns the mean and variance of the Pascal
%   Distribution with parameters N and P.
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      n>=0 (integer)
%      0<=p<=1
%
%   The size of the outputs is the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 2
    error('pascalstat:TooFewInputs',...
        'Requires two input arguments.');
end

% Return NaN for out of range parameters.
n(n<0 | n~=round(n))=NaN;
p(p<0 | p>1)=NaN;

try
    m = n./p;
    v = n.*(1-p)./(p.^2);
catch
    error('pascalstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end