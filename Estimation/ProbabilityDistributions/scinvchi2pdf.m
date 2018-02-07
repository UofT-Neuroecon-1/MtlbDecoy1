function y = scinvchi2pdf(x,v,s)
%SCINVCHI2PDF Scaled Inverse Chi Square probability density function
%   Y = SCINVCHI2PDF(x,v,s) returns the probability density function of 
%   the Scaled Inverse Chi Square Distribution with degrees of freedom V, 
%   scale parameter S^2, at the values in X.
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%        X, V, S > 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011


if nargin < 3
    error('scinvchi2pdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,x,v,s] = distchck(3,x,v,s);

if errorcode > 0
    error('scinvchi2pdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(v,'single') || isa(s,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN for out of range parameters.
y(v<0|s<0|x<0)=NaN;

k=(v>0&s>0&x>0);
if any(k)
    xk=x(k); vk=v(k); sk=s(k);
    term1=(sk.^2).*vk/2;
    term2=(-term1./xk)-gammaln(vk./2);
    term3=(vk./2).*(log(term1)-log(xk))-log(xk);
    y(k)=exp(term2+term3);
end

%round off
y(y<0)=0;


end