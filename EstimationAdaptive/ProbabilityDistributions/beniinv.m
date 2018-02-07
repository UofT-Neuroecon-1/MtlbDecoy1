function x = beniinv(p,a,b,s)
%BENIINV Inverse of the Benini cumulative distribution function
%   X = BENIINV(P,A,B,S) returns the inverse cumulative distribution 
%   function of the Benini distribution with shape parameters A and B 
%   and scale parameter S, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Distribution: Continuous, semi-bounded, [S,Inf)
%   Restrictions:
%        A, B, S > 0
%
%   NOTE: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENIPDF, BENICDF, BENISTAT, BENIFIT, 
%            BENILIKE, BENIRND, BENISF, BENIHAZ
%

%   Mike Sheppard
%   Last Modified 16-May-2012


if nargin ~= 4
    error('beniinv:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode p a b s] = distchck(4,p,a,b,s);

if errorcode > 0
    error('beniinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize x to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single') || isa(s,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<s & s<Inf);
okvar = (0 < p & p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=s(okparam & p==0);
x(okparam & p==1)=Inf;

if any(ok)
    p = p(ok); a = a(ok); b = b(ok); s = s(ok);
    %Solve via transformation
    r=log(1-p);
    q=(-a+sqrt(a.^2-4.*b.*r))./(2.*b);
    x(ok)=s.*exp(q);
end


end
