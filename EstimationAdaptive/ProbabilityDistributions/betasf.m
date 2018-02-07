function s = betasf(x,a,b)
%BETASF Beta survival function.
%   S = BETASF(X,A,B) returns the survival function of the beta 
%   distribution with parameters A and B at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Type: Continuous, bounded, [0,1]
%
%   See also BETAPDF, BETACDF, BETAINV, BETASTAT, 
%            BETAFIT, BETALIKE, BETARND, BETAHAZ
%

%   Mike Sheppard
%   Last Modified 25-Dec-2011


if nargin < 3
    error('betasf:TooFewInputs',...
        'Requires at least three input arguments.');
end

try
    s = 1 - betacdf(x,a,b);
catch
   error('betasf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

