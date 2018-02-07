function [m,v] = asinstat(a,b)
%ASINSTAT Mean and variance for the arcsine distribution.
%   [M,V] = ASINSTAT(A,B) returns the mean and variance of the Arcsine
%   Distribution on the interval [A,B].
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINPDF, ASINCDF, ASININV, ASINFIT, 
%            ASINLIKE, ASINRND, ASINSF, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin == 0
    a = 0;
    b = 1;
end

if nargin == 1,
    error('asinstat:TooFewInputs','Requires either zero or two input arguments.');
end


try
    m = (a+b)/2;
    v = (1/8)*(b-a).^2;
catch err
    error('asinstat:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b);
m(~okparam)=NaN;
v(~okparam)=NaN;

end