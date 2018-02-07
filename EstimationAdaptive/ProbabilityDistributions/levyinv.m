function x = levyinv(p,mu,sigma)
%LEVYINV Inverse of the Lévy cumulative distribution function
%   X = LEVYINV(P,MU,SIGMA) returns the inverse cumulative distribution 
%   function of the Lévy distribution with location parameter MU 
%   and dispersion parameter SIGMA. 
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%     SIGMA > 0
%
%   NOTE: The Lévy Distribution is also known as a van der Waals profile,
%   and if MU=0 it is a special case of the Inverse-Gamma Distribution.
%
%   See also LEVYPDF, LEVYCDF, LEVYSTAT, LEVYFIT, 
%            LEVYLIKE, LEVYRND, LEVYSF, LEVYHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin < 3
   error('levyinv:TooFewInputs','Requires three input arguments.');
end

[errorcode, p,mu,sigma] = distchck(3,p,mu,sigma);

if errorcode > 0
    error('levyinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize p to zero.
if isa(p,'single') || isa(mu,'single') || isa(sigma,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<mu & mu<Inf) & (0<sigma & sigma<Inf);
okvar = (p>0 & p<1);
ok=(okparam & okvar);
x(~okparam)=NaN;
x(okparam & p==0)=mu(okparam & p==0);
x(okparam & p==1)=Inf;

if any(ok)
    p=p(ok); mu=mu(ok); sigma=sigma(ok);
    x(ok)=mu+((sigma./(erfcinv(p).^2))./2);
end



end