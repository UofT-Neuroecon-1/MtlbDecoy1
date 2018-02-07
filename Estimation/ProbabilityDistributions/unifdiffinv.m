function x = unifdiffinv(p)
%UNIFDIFFINV Inverse of the cumulative uniform difference distribution
%   X = UNIFDIFFINV(P) returns the inverse of the uniform difference
%   cumulative distribution at the values in P
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('unifdiffinv:TooFewInputs',...
          'Requires one input argument.'); 
end

%Initialize x to NaN.
if isa(p,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

k0=(p>0)&(p<=(1/2));
k1=(p>1/2)&(p<1);
if any(k0)
    x(k0)=-1+sqrt(2.*p(k0));
end
if any(k1)
    x(k1)=1-sqrt(2.*(1-p(k1)));
end

%Extreme points
x(p==0)=-1;
x(p==1)=1;


end