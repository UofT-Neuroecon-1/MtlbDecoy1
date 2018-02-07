function s = betaprsf(x,a,b)
%BETAPRSF Beta Prime survival function.
%   S = BETAPRSF(X,A,B) returns the survival function of the Beta Prime
%   Distribution with shape parameters A and B, at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRPDF, BETAPRCDF, BETAPRINV, BETAPRSTAT, 
%            BETAPRFIT, BETAPRLIKE, BETAPRRND, BETAPRHAZ,
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin<3,
    error('betaprsf:TooFewInputs','Requires three input arguments.');
end


try
    s=1-betaprcdf(x,a,b);
catch
    error('betaprsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end