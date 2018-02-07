function [m,v]=sknormstat(alpha,mu,sigma)
%SKNORMSTAT Mean and variance of the Skew Normal distribution
%   [M,V]=SKNORMSTAT(ALPHA,MU,SIGMA) returns the mean and variance of the
%   Skew Normal distribution with  shape parameter ALPHA, location
%   parameter MU, and scale parameter SIGMA.
%
%   If SIGMA is not given, then by default SIGMA=1
%   If MU and SIGMA is not given, then by default MU=0, SIGMA=1
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.


%   Mike Sheppard
%   Last Modified 5-Jun-2011

if nargin < 1
    error('sknormstat:TooFewInputs',...
          'Requires at least one input argument.'); 
end

if nargin==1
    mu=0; sigma=1;
end

if nargin==2
    sigma=1;
end

[errorcode alpha mu sigma] = distchck(3,alpha,mu,sigma);

if errorcode > 0
    error('sknormstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(alpha,'single') || isa(mu,'single') || isa(sigma,'single')
    m=zeros(size(mu),'single');
else
    m=zeros(size(mu));
end
v=m;

m=alpha.*mu.*sqrt(2./(pi*(alpha.^2)));
v=(sigma.^2).*(1-((2.*alpha.^2)./(pi.*(1+alpha.^2))));


end