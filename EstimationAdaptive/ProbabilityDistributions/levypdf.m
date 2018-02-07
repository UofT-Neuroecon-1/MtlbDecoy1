function y = levypdf(x,mu,sigma)
%LEVYPDF Lévy probability density function
%   Y = LEVYPDF(X,MU,SIGMA) returns the probability density function of the
%   Lévy Distribution with location parameter MU and dispersion
%   parameter SIGMA.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%      SIGMA > 0
%
%   NOTE: The Lévy Distribution is also known as a van der Waals profile,
%   and if MU=0 it is a special case of the Inverse-Gamma Distribution.
%
%   See also LEVYCDF, LEVYINV, LEVYSTAT, LEVYFIT,
%            LEVYLIKE, LEVYRND, LEVYSF, LEVYHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin ~= 3
   error('levypdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('levypdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<mu & mu<Inf) & (0<sigma & sigma<Inf);
okvar = (x>=mu & x<Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;


if any(ok)
    x=x(ok); sigma=sigma(ok); mu=mu(ok);
    
    term1=sqrt(sigma./(2*pi));
    num=exp(-sigma./(2.*(x-mu)));
    den=(x-mu).^(3/2);
    y(ok) = term1.*(num./den);
end

%Round-off
y(y<0)=0;

end
