function nlogL = baldnichlike(params,data)
%BALDNICHLIKE Negative Balding–Nichols log-likelihood function
%   NLOGL = BALDNICHLIKE(PARAMS,DATA) returns the negative of the
%   Balding-Nichols log-likelihood function for the parameters PARAMS(1)=F
%   and PARAMS(2)=p, given DATA
%
%   Type: Continuous, bounded
%   Restrictions:
%        0 <= F, P <= 1
%
%   See also BALDNICHPDF, BALDNICHCDF, BALDNICHINV, BALDNICHSTAT, 
%            BALDNICHFIT, BALDNICHRND
%

%   Mike Sheppard
%   Last Modified 7-May-2011

if nargin < 2
    error('baldnichlike:TooFewInputs','Requires two input arguments.');
elseif ~isvector(data)
    error('baldnichlike:VectorRequired','DATA must be a vector.');
end
F = params(1);
p = params(2);

%Parametrization of Beta
params2=[(1-F).*p./F,(1-F).*(1-p)./F];
nlogL=betalike(params2,data);
end
