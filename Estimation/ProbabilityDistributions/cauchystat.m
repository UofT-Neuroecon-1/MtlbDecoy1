function [m,v]=cauchystat(a,b)
%CAUCHYSTAT Mean and variance for the Cauchy distribution
%   [M,V] = cauchystat(a,b) returns the mean and variance of the Cauchy
%   Distribution with location parameter A and scale parameter B.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Type: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        B > 0
%
%   Note: The Cauchy Distribution is also known as the Lorentz Distribution
%   Also, the mean and variance are undefined, so function trivially
%   returns NaNs in the proper size.
%
%   See also CAUCHYPDF, CAUCHYCDF, CAUCHYINV, CAUCHYFIT,
%            CAUCHYLIKE, CAUCHYRND, CAUCHYSF, CAUCHYHAZ
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011



if nargin<1, a=0; end
if nargin<2, b=1; end

%No error checking, other then matching sizes, as all outputs are NaNs

try
    m=NaN(size(a+b)); %Common size
    v=m;
catch
    error('cauchystat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end