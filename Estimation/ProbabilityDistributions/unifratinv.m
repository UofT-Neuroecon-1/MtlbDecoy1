function x = unifratinv(p)
%UNIFRATINV Inverse of the Uniform Ratio cumulative distribution
%   Y = UNIFRATINV(P) returns the inverse of the Uniform Ratio cumulative 
%    density at the values in P
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('unifratinv:TooFewInputs',...
          'Requires one input argument.'); 
end

%Initialize x to NaN.
if isa(p,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

k0=(p>=0)&(p<=(1/2));
k1=(p>1/2);
if any(k0)
    x(k0)=2.*p(k0);
end
if any(k1)
    x(k1)=1./(2.*(1-p(k1)));
end


end