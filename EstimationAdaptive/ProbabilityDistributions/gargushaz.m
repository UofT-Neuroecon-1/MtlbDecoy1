function h = gargushaz(x,p,chi,c)
%GARGUSHAZ Generalized ARGUS hazard function
%   H = GARGUSHAZ(X,P,CHI,C) returns the hazard function of the 
%   Generalized ARGUS distribution with power P, curvature CHI,
%   and cut-off C, evaluated at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Type: Continuous, bounded, (0,C)
%   Restrictions:
%        P > -1
%        CHI, C > 0
%
%   Note: The Generalized ARGUS Distribution reduces to the ARGUS
%   distribution with p=1/2.
%
%   See also GARGUSPDF, GARGUSCDF, GARGUSINV, GARGUSSTAT, 
%            GARGUSFIT, GARGUSLIKE, GARGUSRND, GARGUSSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 3
    error('gargushaz:TooFewInputs',...
        'Requires at least three input arguments.');
end

if nargin==3
    c=1;
end

try
    h = garguspdf(x,p,chi,c) ./ gargussf(x,p,chi,c);
catch
    error('gargushaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end