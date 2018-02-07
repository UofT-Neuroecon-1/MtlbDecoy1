function p = cauchycdf(x,a,b)
%CAUCHYPDF Cauchy cumulative distribution function
%   P = CAUCHYCDF(X,A,B) returns the cumulative distribution function
%   of the Cauchy Distribution with location parameter A and scale
%   parameter B, at the values in X. 
%
%   The size of P is the common size of the input arguments. A scalar input
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
%   See also CAUCHYPDF, CAUCHYINV, CAUCHYSTAT, CAUCHYFIT, 
%            CAUCHYLIKE, CAUCHYRND, CAUCHYSF, CAUCHYHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011

if nargin~=1
    error('cauchycdf:TooFewInputs',...
          'Requires at least one input argument.'); 
end    
if nargin==1
    a=0; b=1;
elseif nargin==2
    b=1;
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('cauchycdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to NaN.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (0<b & b<Inf);
okvar = (-Inf < x) & (x < Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x==-Inf)=0;
p(okparam & x==Inf)=1;


if any(ok)
    x = x(ok); a = a(ok); b = b(ok);
    p(ok) = (1/2)+(1/pi)*atan((-a+x)./b);
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end