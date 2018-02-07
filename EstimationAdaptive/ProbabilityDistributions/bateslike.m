function nlogL = bateslike(params,data)
%BATESLIKE Negative log-likelihood for the Bates distribution.
%   NLOGL = BATESLIKE(PARAMS,DATA) returns the negative of the
%   log-likelihood for the Bates distribution, evaluated at parameters
%   PARAMS(1) = N, PARAMS(2) = A, PARAMS(3) = B, given DATA. 
%
%   Type: Continuous, bounded
%   Restrictions:
%        A <= X <= B
%        N > 1        (integer)
%
%   Note: Loss of accuracy occurs when N is greater than 25.
%
%   See also BATESPDF, BATESINV, BATESSTAT, BATESFIT, BATESLIKE, BATESRND
%

%   Mike Sheppard
%   Last Modified: 3-Dec-2011

if nargin < 2
    error('bateslike:TooFewInputs','Requires two input arguments.');
elseif ~isvector(data)
    error('bateslike:VectorRequired','DATA must be a vector.');
end

n = params(1);
a = params(2);
b = params(3);

%Until more accurate method, use definition of log-likelihood
nlogL = -sum(log(batespdf(data(:),n,a,b)));

%Use only real
nlogL = real(nlogL);

end
