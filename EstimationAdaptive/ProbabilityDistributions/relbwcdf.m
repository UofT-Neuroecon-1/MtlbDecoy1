function p = relbwcdf(E,M,G)
%RELBWCDF Relativistic Breit–Wigner cumulative distribution function
%   P = RELBWCDF(E,M,G) returns the cumulative distribution function 
%   of the Relativistic Breit–Wigner distribution with mass of the 
%   resonance M, resonance width (or decay width) G, at the 
%   center-of-mass energy E.
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
%   E,M,G>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 14-Jun-2011


if nargin < 3
    error('relbwcdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,E,M,G] = distchck(3,E,M,G);

if errorcode > 0
    error('relbwcdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(E,'single') || isa(M,'single') || isa(G,'single')
    p=zeros(size(E),'single');
else
    p=zeros(size(E));
end

p(E<0 | M<0 | G<0 )=NaN;


k=(E>=0 & M>=0 & G>=0);
if any(k)
    Ek=E(k); Mk=M(k); Gk=G(k);
    
    %Use transformation of Cauchy CDF
    x=Ek.^2; a=Mk.^2; b=Mk.*Gk;
    catzero=cauchycdf(0,a,b);
    p(k)=(cauchycdf(x,a,b)-catzero)./(1-catzero);
% 
%     
%     %f(E) ~ 1 / [(E^2-M^2)^2 + M^2 G^2)
%     fE=1 ./ (((Ek.^2-Mk.^2).^2) + ((Mk.^2).*(Gk.^2)));
%     %constant of proportionality
%     %Equal to 1 / Integral
%     k1=sqrt(Mk.^2.*(Mk.^2+Gk.^2));
%     num=2*sqrt(2).*Mk.*Gk.*k1;
%     den=sqrt(Mk.^2+k1)*pi;
%     const=num./den;
%     %Normalize so integral is 1
%     y(k)=const.*fE;
end

% %round off
% y(y<0)=0;


end