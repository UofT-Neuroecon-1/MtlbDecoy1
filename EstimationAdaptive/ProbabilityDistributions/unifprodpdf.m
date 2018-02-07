function y = unifprodpdf(x,n)
%UNIFPRODPDF Uniform Product Distribution probability density function
%   Y = UNIFPRODPDF(X,N) returns the probability distribution function of
%   the Uniform Product Distribution of N uniform distributions, evaluated
%   at the values in X
%
%   UNIFPRODPDF(X) returns the Uniform Distribution (N=1)
%
%   Type: Continuous, bounded
%   Restrictions:
%        0 <= X <= 1
%        N>=1         (integer)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 31-Jul-2011


if nargin < 1
    error('unifprodpdf:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin==1, n=1; end


% Return NaN for out of range parameters.
x(x<0 | x>1)=NaN; n(n~=round(n) | n<=0)=NaN;

try
    term1=((-1).^(n-1))./gamma(n);
    term2=(log(x)).^(n-1);
    y=term1.*term2;
catch
    error('unifprodpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end