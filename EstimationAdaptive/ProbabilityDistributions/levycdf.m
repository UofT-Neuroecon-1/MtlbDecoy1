function p = levycdf(x,mu,sigma)
%LEVYCDF Lévy cumulative distribution function
%   P = LEVYCDF(X,MU,SIGMA) returns the cumulative distribution function
%   of the Lévy Distribution with location parameter MU and dispersion
%   parameter SIGMA.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%      SIGMA > 0
%
%   NOTE: The Lévy Distribution is also known as a van der Waals profile,
%   and if MU=0 it is a special case of the Inverse-Gamma Distribution.
%
%   See also LEVYPDF, LEVYINV, LEVYSTAT, LEVYFIT,
%            LEVYLIKE, LEVYRND, LEVYSF, LEVYHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012



if nargin ~= 3
   error('levycdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('levycdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize p to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<mu & mu<Inf) & (0<sigma & sigma<Inf);
okvar = (x>=mu & x<Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<mu)=0;
p(okparam & x==Inf)=1;

if any(ok)
    x=x(ok); sigma=sigma(ok); mu=mu(ok);
    p(ok)=erfc(sqrt(sigma./(2.*(x-mu))));
end

%Catch round off
p(p<0)=0; p(p>1)=1;


end