function y = benipdf(x,a,b,s)
%BENIPDF Benini probability density function
%   Y = BENIPDF(X,A,B,S) returns the probability density function of the 
%   Benini distribution with shape parameters A and B and scale 
%   parameter S, at the values in X.
%
%   The size of Y is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, semi-bounded, [S,Inf)
%   Restrictions:
%        A, B, S > 0
%
%   NOTE: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENICDF, BENIINV, BENISTAT, BENIFIT, 
%            BENILIKE, BENIRND, BENISF, BENIHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012

if nargin ~= 4
    error('benipdf:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode x a b s] = distchck(4,x,a,b,s);

if errorcode > 0
    error('benipdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(s,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<s & s<Inf);
okvar = (s <= x & x < Inf);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x = x(ok); a = a(ok); b = b(ok); s = s(ok);
    y(ok)=((a+2.*b.*log(x./s)).*exp(-a.*log(x./s)-b.*(log(x./s)).^2))./x;
end

%Catch round off
y(y<0)=0;

end
