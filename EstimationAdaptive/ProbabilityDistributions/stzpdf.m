function y = stzpdf(z,n)
%STZPDF Student's Z probability density function
%   Y = STZPDF(z,n) returns the probability density of the
%   Student's Z distribution, with N degrees of freedom, at
%   the values in z.
%
%   Type: Continuous, unbounded
%   Restrictions:
%        N>=1  (integer)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%

%   Mike Sheppard
%   Last Modified 13-Dec-2011


if nargin < 2
    error('stzpdf:TooFewInputs',...
        'Requires at least two input arguments.');
end


[errorcode, z, n] = distchck(2,z,n);


sqrtpi=1.772453850905516;

% Return NaN for out of range parameters.
n(n<1 | n~=round(n))=NaN;

try
    term1=exp(gammaln(n./2)-gammaln((n-1)/2))/sqrtpi;
    term2=(1+(z.^2)).^(-n./2);
    y=term1.*term2;
catch
    error('stzpdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');    
end



end