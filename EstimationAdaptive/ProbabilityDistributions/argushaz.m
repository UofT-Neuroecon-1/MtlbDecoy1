function h = argushaz(x,chi,c)
%ARGUSHAZ ARGUS hazard function
%   H = ARGUSHAZ(X,CHI,C) returns the hazard function of the ARGUS
%   distribution with curvature parameter CHI, cut-off parameter C, 
%   evaluated at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Type: Continuous, bounded, (0,C)
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSINV, ARGUSSTAT, ARGUSFIT, ARGUSLIKE,
%            ARGUSRND, ARGUSSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 2
    error('argushaz:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin ==2
    c=1;
end


try
    h = arguspdf(x,chi,c) ./ argussf(x,chi,c);
catch
   error('argushaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
