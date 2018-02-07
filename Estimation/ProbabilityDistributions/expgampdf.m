function y = expgampdf(x,k,t,mu)
%EXPGAMPDF Exp-Gamma probability density function
%   Y = EXPGAMPDF(X,K,T,MU) returns the probability density function of the
%   Exp-Gamma Distribution with shape parameter K, scale parameter T, 
%   and location parameter MU, at the values of X.
%
%   The size of Y is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        K, T > 0
%
%   Note: The Exp-Gamma Distribution is also known as the Generalized
%   Extreme Value Distribution, and also should not be confused with the
%   Log-Gamma Distribution.
%
%   See also EXPGAMCDF, EXPGAMINV, EXPGAMSTAT, EXPGAMFIT, 
%            EXPGAMLIKE, EXPGAMRND, EXPGAMSF, EXPGAMHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 4
    error('expgampdf:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode x k t mu] = distchck(4,x,k,t,mu);

if errorcode > 0
    error('expgampdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to NaN.
if isa(x,'single') || isa(k,'single') || isa(t,'single') || isa(mu,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<k & k<Inf) & (0<t & t<Inf) & (-Inf<mu & mu<Inf);
okvar = (-Inf<x & x<Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    x = x(ok); k = k(ok); t = t(ok); mu = mu(ok);
    %For accuracy use gammaln
    lognum=(-exp((-mu+x)./t))+(k.*(-mu+x)./t);
    logden=gammaln(k)+log(t);
    y(ok)=exp(lognum-logden);
end

%Catch round off
y(y<0)=0;

end