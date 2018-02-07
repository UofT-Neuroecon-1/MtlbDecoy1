function y = cauchypdf(x,a,b)
%CAUCHYPDF Cauchy probability density function
%   Y = CAUCHYPDF(X,A,B) returns the probability density function of the
%   Cauchy Distribution with location parameter A and scale parameter B,
%   at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Distribution: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        B > 0
%
%   Note: The Cauchy Distribution is also known as the Lorentz Distribution
%
%   See also CAUCHYCDF, CAUCHYINV, CAUCHYSTAT, CAUCHYFIT,
%            CAUCHYLIKE, CAUCHYRND, CAUCHYSF, CAUCHYHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011

if nargin ~= 1
    error('cauchypdf:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin==1
    a=0; b=1;
elseif nargin==2
    b=1;
end


[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('cauchypdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (0<b & b<Inf);
okvar = (-Inf < x) & (x < Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok),
    x=x(ok); a=a(ok); b=b(ok);
    y(ok)=(b.*pi.*(1+((-a+x)./b).^2)).^(-1);
end

%Catch round off
y(y<0)=0;

end