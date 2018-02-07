function [m,v]=nigstat(mu,alpha,beta,delta)
%NIGSTAT Mean and variance for the Normal-Inverse Gaussian (NIG)
%   [M,V] = NIGSTAT(MU,ALPHA,BETA,DELTA) returns the mean and variance
%   for the Normal-Inverse Gaussian (NIG) distirbution with location MU,
%   tail heavyness ALPHA, asymmetry parameter BETA, and scale parameter
%   DELTA.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%   alpha>beta
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs. 
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 4
    error('nigstat:TooFewInputs',...
          'Requires four input argument.');
end


[errorcode, mu,alpha,beta,delta] = distchck(4,mu,alpha,beta,delta);

if errorcode > 0
    error('nigstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(mu,'single') || isa(alpha,'single') || isa(beta,'single') || isa(delta,'single')
    m=zeros(size(mu),'single');
else
    m=zeros(size(mu));
end
v=m;

gamma=sqrt(alpha.^2-beta.^2);
m=mu+(delta.*beta)./gamma;
v=delta.*(alpha.^2)./(gamma.^3);
end

