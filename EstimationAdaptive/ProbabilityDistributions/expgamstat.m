function [m,v]=expgamstat(k,t,mu)
%EXPGAMSTAT Mean and variance for the Exp-Gamma distribution
%   [M,V] = EXPGAMSTAT(K,T,MU) returns the mean and variance of the
%   Exp-Gamma Distribution with shape parameter K, scale parameter
%   T, and location parameter MU.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        K, T > 0
%
%   Note: The Exp-Gamma Distribution is also known as the Generalized
%   Extreme Value Distribution, and also should not be confused with the
%   Log-Gamma Distribution.
%
%   See also EXPGAMPDF, EXPGAMCDF, EXPGAMINV, EXPGAMFIT, 
%            EXPGAMLIKE, EXPGAMRND, EXPGAMSF, EXPGAMHAZ
%

%   Mike Sheppard
%   Last Modified 10-May-2011

if nargin < 3
    error('expgamstat:TooFewInputs',...
        'Requires at least four input argument.');
end

[errorcode k t mu] = distchck(3,k,t,mu);

if errorcode > 0
    error('expgamstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize M and V to zero.
if isa(k,'single') || isa(t,'single') || isa(mu,'single')
    m=zeros(size(k),'single');
else
    m=zeros(size(k));
end
v=m;

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<k & k<Inf) & (0<t & t<Inf) & (-Inf<mu & mu<Inf);
ok=(okparam);
m(~ok)=NaN; v(~ok)=NaN;

if any(ok)
    k=k(ok); t=t(ok); mu=mu(ok);
    m(ok)=mu+t.*psi(0,k);
    v(ok)=(t.^2).*psi(1,k);
end

end