function [m,v] = relbwstat(M,G)
%RELBWSTAT Mean and variance for the Relativistic Breit–Wigner distribution
%   [M,V] = RELBWSTAT(M,G) returns the mean and variance of the Relativistic
%   Breit–Wigner probability density function with mass of the resonance M,
%   resonance width (or decay width) G.
%
%   This uses natural units of h_bar=c=1
%
%   f(E) ~ 1 / [(E^2-M^2)^2 + M^2 G^2]
%
%   Note: The (non-relativistic) Breit–Wigner distribution is the Cauchy
%   Distribution
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%   M,G>0
%
%   The size of outputs  is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 2
    error('relbwstat:TooFewInputs',...
        'Requires two input arguments.');
end

try
    %Expand size if necessary
    M=M+zeros(size(G));
    G=G+zeros(size(M));
catch
    error('relbwstat:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end



% Initialize Y to zero.
if isa(M,'single') || isa(G,'single')
    m=zeros(size(M),'single');
else
    m=zeros(size(M));
end
v=m;

m(M<0 | G<0 )=NaN;
v(M<0 | G<0 )=NaN;


k=(M>=0 & G>=0);
if any(k)
    Mk=M(k); Gk=G(k);
    %Expected value of integral, without constant of proportionality
    term1=(pi+2.*atan(Mk./Gk))./(4.*Gk.*Mk);
    %constant of proportionality
    %Equal to 1 / Integral
    k1=sqrt(Mk.^2.*(Mk.^2+Gk.^2));
    num=2*sqrt(2).*Mk.*Gk.*k1;
    den=sqrt(Mk.^2+k1)*pi;
    const=num./den;
    %Normalize so integral is 1
    m(k)=const.*term1;
end


v=NaN(size(m));  %Unknown formula for variance

end
