function s = argussf(x,chi,c)
%ARGUSSF ARGUS survival function
%   S = ARGUSSF(X,CHI,C) returns the survival function of the ARGUS
%   distribution with curvature parameter CHI, cut-off parameter C, 
%   evaluated at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Type: Continuous, bounded, (0,C)
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSINV, ARGUSSTAT, 
%            ARGUSFIT, ARGUSLIKE, ARGUSRND, ARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 2
    error('argussf:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin ==2
    c=1;
end


try
    s=1-arguscdf(x,chi,c);
catch
   error('argussf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
