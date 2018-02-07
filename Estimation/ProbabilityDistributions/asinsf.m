function s = asinsf(x,a,b)
%ASINSF Arcsine survival function
%   S = ASINSF(X,A,B) returns the survival function of the
%   Arcsine Distribution on the interval [A,B] at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINPDF, ASINCDF, ASININV, ASINSTAT, 
%            ASINFIT, ASINLIKE, ASINRND, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if (nargin<1)
    error('asinsf:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin == 1
    a=0; b=1;
elseif nargin==2
    b=1;
end


try
    s=1-asincdf(x,a,b);
catch err
   error('asinsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end