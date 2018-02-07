function x = invgaminv(p,a,b)
%INVGAMINV Inverse of the Inverse-Gamma cumulative distribution function
%   X = INVGAMINV(P,A,B) returns the inverse cumulative distribution 
%   function of the Inverse-Gamma distribution with shape parameter A 
%   and scale parameter B, evaluated at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A, B > 0
%
%   See also INVGAMPDF, INVGAMCDF, INVGAMSTAT, INVGAMFIT, 
%            INVGAMLIKE, INVGAMRND, INVGAMSF, INVGAMHAZ
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 3
   error('invgaminv:TooFewInputs','Requires three input arguments.');
end

[errorcode, p,a,b] = distchck(3,p,a,b);

if errorcode > 0
    error('invgaminv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end


k=(a>0 & b>0 & p>0 & p<1);
if any(k)
    x(k)=b(k)./gammaincinv(p(k),a(k),'upper');
end

%return zero or NaN for out of range or invalid
x(a>0 & b>0 & p==0)=0;
x(a>0 & b>0 & p==Inf)=Inf;
x(a<0 | b<0 | p<0 | p>1) = NaN;

end