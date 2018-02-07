function y = unifratpdf(x)
%UNIFRATPDF Uniform Ratio Distribution probability density function
%   Y = UNIFRATPDF(X) returns the probability density function of the
%   Uniform Ratio Distribution at the values in X
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('unifratpdf:TooFewInputs',...
          'Requires one input argument.'); 
end

%Initialize y to zero.
if isa(x,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

k0=(x>0)&(x<=1); k1=(x>1);
if any(k0)
    y(k0)=1/2;
end
if any(k1)
    y(k1)=1./(2.*(x(k1).^2));
end

y(x<0)=0;

end