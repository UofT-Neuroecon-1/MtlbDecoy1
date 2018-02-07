function h = gbetaprhaz(x,a,b,p,q)
%GBETAHAZ Generalized Beta Prime hazard function.
%   H = GBETAPRHAZ(X,A,B,P,Q) returns the hazard function of the 
%   Generalized Beta Prime Distribution with shape parameters
%   A, B, and P, and scale parameter Q at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B, P, Q > 0
%
%   Note: The Generalized Beta Prime Distribution is also known as the
%   Generalized Inverted Beta Distribution, or Generalized Beta
%   Distribution of the Second Kind
%
%   See also GBETAPRPDF, GBETAPRCDF, GBETAPRINV, GBETAPRSTAT, 
%            GBETAPRFIT, GBETAPRLIKE, GBETAPRRND, GBETAPRSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 5
    error('gbetaprhaz:TooFewInputs','Requires five input arguments.');
end

try
    h = gbetaprpdf(x,a,b,p,q) ./ gbetaprsf(x,a,b,p,q);
catch
    error('gbetaprhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end