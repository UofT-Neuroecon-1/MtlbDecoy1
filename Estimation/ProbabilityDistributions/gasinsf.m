function s = gasinsf(x,alpha,a,b)
%GASINSF Generalized arcsine survival function.
%   S = GASINSF(X,ALPHA,A,B) returns the survival function of the 
%   Generalized Arcsine Distribution with shape parameter ALPHA
%   on the interval [A,B] at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
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
%            GASINRND, GASINHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if (nargin~=2)&&(nargin~=4)
    error('gasinsf:TooFewInputs',...
        'Requires either two or four input arguments.');
end

if nargin == 2
    a = 0;
    b = 1;
end

try
    s = 1-gasincdf(x,alpha,a,b);
catch
    error('gasinsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end