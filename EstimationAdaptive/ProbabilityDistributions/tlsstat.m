function [m,v] = tlsstat(mu,sigma,nu)
%TLSSTAT Mean and variance for the T Location-Scale Distribution.
%   [M,V] = TLSSTAT(MU,SIGMA,NU) returns the mean and variance of the
%   T Location-Scale Distribution with location parameter MU, scale 
%   parameter SIGMA, and NU degrees of freedom.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA, NU > 0
%
%   See also TLSPDF, TLSCDF, TLSINV, TLSFIT, TLSLIKE, TLSRND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 3
    error('tlsstat:TooFewInputs',...
        'Requires three input arguments.');
end


[errorcode,mu,sigma,nu] = distchck(3,mu,sigma,nu);

if errorcode > 0
    error('tlsstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%pre-allocate memory
m=zeros(size(mu)); v=m;

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;


%Modified from addtls.m in statistics toolbox
%-----

k=(nu<=1);
m(k)=NaN;
m(~k)=mu(~k);

k=(nu<=2);
v(k)=Inf;
v(~k)=sigma(~k).^2 .* nu(~k) ./ (nu(~k) - 2);

%-----


end