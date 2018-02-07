function x = unifprodinv(p,n)
%UNIFPRODINV Cumulative Uniform Product probability density
%   X = UNIFPRODINV(P,N) returns the Inverse Cumulative Uniform Product 
%   probability density of N uniform distributions, at the values in P
%
%   Default value for n is 1.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 13-Dec-2011


if nargin < 1
    error('unifprodinv:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin==1, n=1; end


% Return NaN for out of range parameters.
p(p<0 | p>1)=NaN; n(n~=round(n))=NaN;

try
    x=exp(-gammaincinv(p,n,'upper'));
catch
     error('unifprodinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end



end