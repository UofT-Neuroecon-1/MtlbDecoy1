function y = relbwpdf(E,M,G)
%RELBWPDF Relativistic Breit–Wigner probability density function
%   Y = RELBWPDF(E,M,G) returns the probability density function of the
%   Relativistic Breit–Wigner Distribution with mass of the resonance M,
%   resonance width (or decay width) G, at the center-of-mass energy E.
%
%   This uses natural units of h_bar = c = 1
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%        E, M, G > 0
%
%   The Relativistic Breit-Wigner Distribution is proportional to
%   f(E) ~ 1 / [(E^2-M^2)^2 + M^2 G^2]
%
%   Note: The (non-relativistic) Breit–Wigner distribution is the Cauchy
%   Distribution
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011


if nargin < 3
    error('relbwpdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,E,M,G] = distchck(3,E,M,G);

if errorcode > 0
    error('relbwpdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(E,'single') || isa(M,'single') || isa(G,'single')
    y=zeros(size(E),'single');
else
    y=zeros(size(E));
end

y(E<0 | M<0 | G<0 )=NaN;


k=(E>=0 & M>=0 & G>=0);
if any(k)
    Ek=E(k); Mk=M(k); Gk=G(k);
    %f(E) ~ 1 / [(E^2-M^2)^2 + M^2 G^2)
    fE=1 ./ (((Ek.^2-Mk.^2).^2) + ((Mk.^2).*(Gk.^2)));
    %constant of proportionality
    %Equal to 1 / Integral
    k1=sqrt(Mk.^2.*(Mk.^2+Gk.^2));
    num=2*sqrt(2).*Mk.*Gk.*k1;
    den=sqrt(Mk.^2+k1)*pi;
    const=num./den;
    %Normalize so integral is 1
    y(k)=const.*fE;
end

%round off
y(y<0)=0;


end