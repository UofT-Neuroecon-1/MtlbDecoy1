function h = betaprhaz(x,a,b)
%BETAPRHAZ Beta Prime hazard function.
%   H = BETAPRHAZ(X,A,B) returns the hazard function of the Beta Prime
%   Distribution with shape parameters A and B, at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRPDF, BETAPRCDF, BETAPRINV, BETAPRSTAT, BETAPRFIT,
%            BETAPRLIKE, BETAPRRND, BETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin<3,
    error('betaprhaz:TooFewInputs','Requires three input arguments.');
end


try
    h = betaprpdf(x,a,b) ./ betaprsf(x,a,b);
catch
    error('betaprhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end