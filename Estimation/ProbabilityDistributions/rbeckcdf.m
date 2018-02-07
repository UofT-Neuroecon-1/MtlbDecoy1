function p = rbeckcdf(x,s1,s2)
%RBECKCDF Reduced Beckmann cumulative distribution function
%   P = RBECKCDF(X,S1,S2) returns the cumulative distribution function of 
%   the Reduced Beckmann distribution with standard deviations S1 and S2, 
%   at the values of X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The full Beckmann Distribution
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
%   See also RBECKPDF, RBECKINV, RBECKSTAT, RBECKFIT, 
%            RBECKLIKE, RBECKRND, RBECKSF, RBECKHAZ
%

%   Mike Sheppard
%   Last Modified 3-Dec-2011


if nargin ~= 3
    error('rbeckcdf:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode x s1 s2] = distchck(3,x,s1,s2);

if errorcode > 0
    error('rbeckcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(s1,'single') || isa(s2,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<s1 & s1<Inf) & (0<s2 & s2<Inf);
%okvar = (0 <= x & x < Inf);
s1(~okparam)=NaN; s2(~okparam)=NaN;

k=(s1==s2);
if any(k)
    %Simplified expression
    p(k) = 1 - exp((-x(k).^2)./(2*s2(k).^2));
end
if any(~k)
    xk=x(~k);
    s1k=s1(~k);
    s2k=s2(~k);
    mk=min(s1k,s2k);
    Mk=max(s1k,s2k);
    comk=(sqrt(1-((mk./Mk).^4)).*Mk.*xk)./(2.*mk.*sqrt(s1k.^2+s2k.^2));
    term1=sqrt(1+(2*mk./(-mk+Mk)));
    term2=sqrt(1-(2*mk./(mk+Mk)));
    %Use MarcumQ2 function, derivation of MarcumQ in signal processing
    %using statistics toolbox
    a1=comk.*term1; b1=comk.*term2;
    a2=comk.*term2; b2=comk.*term1;
    p(~k)=marcumq2(a1,b1)-marcumq2(a2,b2);
end

%Edge cases
p(okparam & x<=0)=0;
p(okparam & x==Inf)=1;

%Catch round off
p(p<0)=0; p(p>1)=1;

end


function Q=marcumq2(a,b)
%MARCUMQ function is given in the Signal Processing Toolbox, but the only
%toolbox we are assuming is the Statistics Toolbox.
%The Marcum Q function with parameters [A,B,M] can be derived from the CDF
%of the non-central chi-square distribution
%Q_m(a,b) = 1 - ncx2cdf(b^2,2*m,a^2)
Q = 1 - ncx2cdf(b.^2,2,a.^2);
end