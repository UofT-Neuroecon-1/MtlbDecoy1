function h = benihaz(x,a,b,s)
%BENIHAZ Benini hazard function
%   H = BENIHAZ(X,A,B,S) returns the Benini hazard function with 
%   shape parameters A and B and scale parameter S, at the values in X.
%
%   The size of H is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, semi-bounded, [S,Inf)
%   Restrictions:
%        A, B, S > 0
%
%   Note: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENIPDF, BENICDF, BENIINV, BENISTAT, BENIFIT, BENILIKE,
%            BENIRND, BENISF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 4
    error('benihaz:TooFewInputs',...
          'Requires at least four input argument.'); 
end

[errorcode x a b s] = distchck(4,x,a,b,s);

if errorcode > 0
    error('benihaz:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize h to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(s,'single')
    h=zeros(size(x),'single');
else
    h=zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<s & s<Inf);
okvar = (s <= x & x < Inf);
ok=(okparam & okvar);
h(~ok)=NaN;

%Could use definition of hazard function, but direct formula is less computational
if any(ok),
    x = x(ok); a = a(ok); b = b(ok); s = s(ok);
    h(ok) = (a+2.*b.*log(x./s))./x;
end


end