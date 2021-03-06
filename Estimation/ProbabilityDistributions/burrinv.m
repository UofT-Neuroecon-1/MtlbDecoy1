function x = burrinv(p,c,k)
%BURRINV Inverse of the Burr cumulative distribution function (cdf)
%   X = BURRINV(P,C,K) returns the inverse cumulative distribution
%   function of the Burr distribution with shape parameters C and K, 
%   at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%      C, K > 0
%
%   NOTE: The Burr Distribution is also known as the Burr Type XII
%   Distribution, and is a special case of the Sing-Maddala Distribution.
%   The notation for the distribution used here is F(X<x) = 1 - (1+x^C)^-K   
%
%   See also BURRPDF, BURRCDF, BURRSTAT, BURRFIT, 
%            BURRLIKE, BURRRND, BURRSF, BURRHAZ
%


%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin ~= 3
    error('burrinv:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode p c k] = distchck(3,p,c,k);

if errorcode > 0
    error('burrinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize x to zero.
if isa(p,'single') || isa(c,'single') || isa(k,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<k & k<Inf);
okvar = (0<p & p<1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=0;
x(okparam & p==1)=Inf;


if any(ok),
    p = p(ok); c = c(ok); k = k(ok);
    x(ok) =(((1-p).^(-1./k))-1).^(1./c);
end



end

