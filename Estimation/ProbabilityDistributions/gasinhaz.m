function h = gasinhaz(x,alpha,a,b)
%GASINHAZ Generalized arcsine hazard function.
%   H = GASINHAZ(X,ALPHA,A,B) returns the hazard function of the 
%   Generalized Arcsine Distribution with shape parameter ALPHA
%   on the interval [A,B] at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Type: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        0 < ALPHA < 1
%
%   See also GASINPDF, GASINCDF, GASININV, GASINSTAT, GASINFIT, GASINLIKE,
%            GASINRND, GASINSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if (nargin~=2)&&(nargin~=4)
    error('gasinhaz:TooFewInputs',...
        'Requires either two or four input arguments.');
end

if nargin == 2
    a = 0;
    b = 1;
end

try
    h = gasinpdf(x,alpha,a,b) ./ gasinsf(x,alpha,a,b);
catch
    error('gasinsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end