function y = invgampdf(x,a,b)
%INVGAMPDF Inverse-Gamma probability density function
%   Y = INVGAMPDF(X,A,B) returns the probability density function of the
%   Inverse-Gamma Distribution with shape parameter A and scale 
%   parameter B, evaluated at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%      A, B > 0
%
%   See also INVGAMCDF, INVGAMINV, INVGAMSTAT, INVGAMFIT,
%            INVGAMLIKE, INVGAMRND, INVGAMSF, INVGAMHAZ
%

%   Mike Sheppard
%   Last Modified 3-Jul-2011


if nargin ~= 3
   error('invgampdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,a,b] = distchck(3,x,a,b);

if errorcode > 0
    error('invgampdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf);
okvar = (0<=x & x<Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    a=a(ok); b=b(ok); x=x(ok);
    lognum=(-b./x)+(a.*(log(b)-log(x)));
    logden=log(x)+gammaln(a);
    y(ok)=exp(lognum-logden);
end

%Round-off
y(y<0)=0;


end