function [m,v] = gbetaprstat(a,b,p,q)
%GBETAPRSTAT Mean and variance for the Generalized Beta Prime Distribution.
%   [M,V] = GBETAPRSTAT(A,B,P,Q) returns the mean and variance of the 
%   Generalized Beta Prime distribution with shape parameters 
%   A, B, and P, and scale parameter Q.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B, P, Q > 0
%
%   Note: The Generalized Beta Prime Distribution is also known as the
%   Generalized Inverted Beta Distribution, or Generalized Beta
%   Distribution of the Second Kind
%
%   See also GBETAPRPDF, GBETAPRCDF, GBETAPRINV, GBETAPRFIT, 
%            GBETAPRLIKE, GBETAPRRND, GBETAPRSF, GBETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin < 4
    error('gbetaprstat:TooFewInputs',...
          'Requires at least four input arguments.'); 
end

[errorcode a b p q] = distchck(4,a,b,p,q);

if errorcode > 0
    error('gbetaprstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


if isa(a,'single') || isa(b,'single') || isa(p,'single') || isa(q,'single')
   m = zeros(size(a),'single');
else
   m = zeros(size(a));
end
v = m;

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<p & p<Inf) & (0<q & q<Inf);

%Use gammaln and exp in place of gamma
m=q.*exp(gammaln(a+(1./p))+gammaln(b-(1./p))-gammaln(a)-gammaln(b));
secondmoment=(q.^2).*exp(gammaln(a+(2./p))+gammaln(b-(2./p))-gammaln(a)-gammaln(b));
v=secondmoment-m.^2;

m(~okparam)=NaN; v(~okparam)=NaN;
m(okparam & b.*p<=1)=Inf;
v(okparam & b.*p<=2)=NaN;


end