function p = unifratcdf(x)
%UNIFRATCDF Cumulative Uniform Ratio probability density
%   P = UNIFRATCDF(X) returns the cumulative Uniform Ratio density at the
%   values in X
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('unifratcdf:TooFewInputs',...
          'Requires one input argument.'); 
end

%Initialize p to NaN.
if isa(x,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

k0=(x>0)&(x<=1);
k1=(x>1);
if any(k0)
    p(k0)=x(k0)/2;
end
if any(k1)
    p(k1)=(-1+(2.*x(k1)))./(2.*x(k1));
end


end