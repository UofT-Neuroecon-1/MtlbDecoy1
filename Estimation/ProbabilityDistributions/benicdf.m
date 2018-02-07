function p = benicdf(x,a,b,s)
%BENICDF Benini cumulative distribution function
%   P = BENICDF(X,A,B,S) returns the cumulative distribution function
%   of the Benini distribution with shape parameters A and B, and scale
%   parameter S, at the values in X.
%
%   The size of P is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, semi-bounded, [S,Inf)
%   Restrictions:
%        A, B, S > 0
%
%   Note: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENIPDF, BENIINV, BENISTAT, BENIFIT, 
%            BENILIKE, BENIRND, BENISF, BENIHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011

if nargin ~= 4
    error('benicdf:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode x a b s] = distchck(4,x,a,b,s);

if errorcode > 0
    error('benicdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(s,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<s & s<Inf);
okvar = (s <= x & x < Inf);
ok = (okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<s)=0;
p(okparam & x==Inf)=1;

if any(ok)
    x = x(ok); a = a(ok); b = b(ok); s = s(ok);
    p(ok)=1-exp(-a.*log(x./s)-b.*(log(x./s)).^2);
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end