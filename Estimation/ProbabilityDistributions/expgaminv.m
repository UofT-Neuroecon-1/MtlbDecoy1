function x = expgaminv(p,k,t,mu)
%EXPGAMINV Inverse of the Exp-Gamma cumulative distribution function
%   X = EXPGAMINV(P,K,T,MU) returns the inverse cumulative distribution 
%   function of the Exp-Gamma Distribution with shape parameter K, 
%   scale parameter T, and location parameter MU at the values in P.
%
%   The size of X is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        K, T > 0
%
%   NOTE: The Exp-Gamma Distribution is also known as the Generalized
%   Extreme Value Distribution, and also should not be confused with the
%   Log-Gamma Distribution.
%
%   See also EXPGAMPDF, EXPGAMCDF, EXPGAMSTAT, EXPGAMFIT, 
%            EXPGAMLIKE, EXPGAMRND, EXPGAMSF, EXPGAMHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 4
    error('expgaminv:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode p k t mu] = distchck(4,p,k,t,mu);

if errorcode > 0
    error('expgaminv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to NaN.
if isa(p,'single') || isa(k,'single') || isa(t,'single') || isa(mu,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<k & k<Inf) & (0<t & t<Inf) & (-Inf<mu & mu<Inf);
okvar = (0<p & p<1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=-Inf;
x(okparam & p==1)=Inf;

if any(ok)
    p = p(ok); k = k(ok); t = t(ok); mu = mu(ok);
    x(ok) = mu+t.*log(gammaincinv(p,k));
end

end