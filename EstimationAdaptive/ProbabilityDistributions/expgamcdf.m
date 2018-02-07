function p = expgamcdf(x,k,t,mu)
%EXPGAMCDF Exp-Gamma cumulative distribution function
%   P = EXPGAMCDF(X,K,T,MU) returns the cumulative distribution function
%   of the Exp-Gamma Distribution with shape parameter K, scale 
%   parameter T, and location parameter MU, at the values of X.
%
%   The size of P is the common sizes of the inputs. A scalar input   
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
%   See also EXPGAMPDF, EXPGAMINV, EXPGAMSTAT, EXPGAMFIT, 
%            EXPGAMLIKE, EXPGAMRND, EXPGAMSF, EXPGAMHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 4
    error('expgamcdf:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode x k t mu] = distchck(4,x,k,t,mu);

if errorcode > 0
    error('expgamcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to NaN.
if isa(x,'single') || isa(k,'single') || isa(t,'single') || isa(mu,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<k & k<Inf) & (0<t & t<Inf) & (-Inf<mu & mu<Inf);
okvar = (-Inf<x & x<Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x==-Inf)=0;
p(okparam & x==Inf)=1;


if any(ok)
    x = x(ok); k = k(ok); t = t(ok); mu = mu(ok);
    term=exp((x-mu)./t);
    %http://reference.wolfram.com/mathematica/ref/ExpGammaDistribution.html
    %CDF = GammaRegularized[kk,0,term] [Mathematica]
    %    = (Gamma[kk,0]-Gamma[kk,term])/Gamma[kk] [Mathematica]
    %Notation: 
    %   Gamma[a,b] [Mathematica] = (1-gammainc(b,a))*gamma(a)  [MATLAB]
    %CDF = (((1-gammainc(0,kk))*gamma(kk))-((1-gammainc(term,kk))*gamma(kk)))/gamma(kk)
    %    = gammainc(term,kk)  [MATLAB]
    p(ok)=gammainc(term,k);
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end