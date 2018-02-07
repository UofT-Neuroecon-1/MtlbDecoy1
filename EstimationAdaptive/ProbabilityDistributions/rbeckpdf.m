function y = rbeckpdf(x,s1,s2)
%RBECKPDF (Reduced) Beckmann probability density function
%   Y = RBECKPDF(X,S1,S2) returns the probability density function of the
%   (reduced) Beckmann Distribution with standard deviations S1 and S2,
%   at the values of X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The full Beckmann Distribution:
%   If {x,y} follows the Bivariate Normal Distribution with means {u1,u2},
%   standard deviations {s1,s2} and correlation rho, then sqrt[x^2+y^2]
%   follows the Beckmann Distribution [u1,u2,s1,s2,rho]
%
%   The reduced Beckmann distribution assumes u1=u2=rho=0 so the bivariate
%   normal distributions are centered at (0,0) with covariance matrix
%   [s1^2 0; 0 s2^2]
%
%   Distribution: Continuous, semi-bounded, [0,Inf) 
%   Restrictions:
%        S1, S2 > 0
%
%   See also RBECKCDF, RBECKINV, RBECKSTAT, RBECKFIT, 
%            RBECKLIKE, RBECKRND, RBECKSF, RBECKHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin ~= 3
    error('rbeckpdf:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode x s1 s2] = distchck(3,x,s1,s2);

if errorcode > 0
    error('rbeckpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(s1,'single') || isa(s2,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<s1 & s1<Inf) & (0<s2 & s2<Inf);
okvar = (0 <= x & x < Inf);
s1(~okparam)=NaN; 
s2(~okparam)=NaN;
k=(s1==s2);

if any(k)
    %Simplified expression
    xk=x(k); s2k=s2(k);
    y(k) = exp((-xk.^2)./(2*s2k.^2)).*xk./(s2k.^2);
end
if any(~k)
    xk=x(~k); s1k=s1(~k); s2k=s2(~k);
    term1=(-((s1k.^2)+(s2k.^2)).*(xk.^2))./(4*(s1k.^2).*(s2k.^2));
    term2=(1/4).*(-(s1k.^-2)+(s2k.^-2)).*(xk.^2);
    y(~k)=exp(term1).*xk.*besseli(0,term2)./(s1k.*s2k);
end

%Edge cases
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

%Catch round off
y(y<0)=0;

end