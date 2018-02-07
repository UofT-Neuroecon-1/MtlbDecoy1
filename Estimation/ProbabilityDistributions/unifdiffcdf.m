function p = unifdiffcdf(x)
%UNIFDIFFCDF Cumulative Uniform Difference Distribution
%   P = UNIFDIFFCDF(X) returns the cumulative uniform difference
%   distribution at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('unifdiffcdf:TooFewInputs',...
          'Requires one input argument.'); 
end

%Initialize p to NaN.
if isa(x,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

%Three ranges
k0=(x>-1)&(x<=0); k1=(x>0)&(x<1); k2=(x>=1);
if any(k0)
    p(k0)=(1/2).*((1+x(k0)).^2);
end
if any(k1)
    p(k1)=(1/2)+x(k1)-((x(k1).^2)/2);
end
p(k2)=1;


end