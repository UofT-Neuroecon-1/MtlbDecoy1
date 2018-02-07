function h = betahaz(x,a,b)
%BETAHAZ Beta survival function.
%   H = BETAHAZ(X,A,B) returns the hazard function of the beta 
%   distribution with parameters A and B at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Type: Continuous, bounded, [0,1]
%
%   See also BETAPDF, BETACDF, BETAINV, BETASTAT, 
%            BETAFIT, BETALIKE, BETARND, BETASF
%

%   Mike Sheppard
%   Last Modified 25-Dec-2011


if nargin < 3
    error('betahaz:TooFewInputs',...
        'Requires at least three input arguments.');
end

try
    h = betapdf(x,a,b) ./ betasf(x,a,b);
catch
   error('betahaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

