function p = invgamcdf(x,a,b)
%INVGAMCDF Inverse-Gamma cumulative distribution function
%   P = INVGAMCDF(X,A,B) returns the cumulative distribution function 
%   of the Inverse-Gamma Distribution with shape parameter A and scale
%   parameter B, evaluated at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%      A, B > 0
%
%   See also INVGAMPDF, INVGAMINV, INVGAMSTAT, INVGAMFIT,
%            INVGAMLIKE, INVGAMRND, INVGAMSF, INVGAMHAZ
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin ~= 3
   error('invgamcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,a,b] = distchck(3,x,a,b);

if errorcode > 0
    error('invgamcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf);
okvar = (0 < x) & (x < Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<=0)=0;
p(okparam & x==Inf)=1;


if any(ok)
    x=x(ok); a=a(ok); b=b(ok);
    p(ok)=gammainc(b./x,a,'upper');
end

%Catch round off
p(p<0)=0; p(p>1)=1;


end