function [m,v]=maxwstat(a)
%MAXWSTAT Mean and variance for the Maxwell-Boltzmann distribution
%   [M,V] = MAXWSTAT(a) returns the mean and variance of the
%   Maxwell-Boltzmann distribution with scale parameter A.

%   Mike Sheppard
%   Last Modified 24-Mar-2011


if nargin < 1
    error('maxwstat:TooFewInputs',...
          'Requires at least one input arguments.'); 
end

% Initialize y to NaN.
if isa(a,'single')
   m = zeros(size(a),'single');
else
   m = zeros(size(a));
end
v=m;

k=(a>0);
if any(k)
m(k) = 2*a(k)*sqrt(2/pi);
v(k) = (a(k)).^2.*(-8+3*pi)./pi;
end


end