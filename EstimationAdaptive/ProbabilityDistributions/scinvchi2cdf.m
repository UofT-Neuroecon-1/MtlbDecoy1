function p = scinvchi2cdf(x,v,s)
%SCINVCHI2CDF Scaled Inverse Chi Square cumulative distribution function
%   P = SCINVCHI2CDF(x,v,s) returns the cumulative distribution function 
%   of the Scaled Inverse Chi Square distribution with degrees of 
%   freedom V, and scale parameter S^2.
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%      x,v,s>0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011


if nargin < 3
    error('scinvchi2cdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,x,v,s] = distchck(3,x,v,s);

if errorcode > 0
    error('scinvchi2cdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize P to zero.
if isa(x,'single') || isa(v,'single') || isa(s,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

p(v<0|s<0|x<0)=NaN;

k=(v>0&s>0&x>0);
if any(k)
    xk=x(k); vk=v(k); sk=s(k);
    term1=(sk.^2).*vk/(2.*x);
    p(k)=gammainc(term1,vk./2,'upper');
end



end