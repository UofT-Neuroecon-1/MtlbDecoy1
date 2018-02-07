function p = unifprodcdf(x,n)
%UNIFPRODCDF Cumulative Uniform Product probability density
%   P = UNIFPRODCDF(X,N) returns the Cumulative Uniform Product 
%   probability density of N uniform distributions, at the values in X
%
%   Default value for n is 1.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 13-Dec-2011


if nargin < 1
    error('unifprodcdf:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin==1, n=1; end

% Return NaN for out of range parameters.
x(x<0 | x>1)=NaN; n(n~=round(n))=NaN;

try
    p=1-gammainc(-log(x),n);
catch
    error('unifprodcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


end