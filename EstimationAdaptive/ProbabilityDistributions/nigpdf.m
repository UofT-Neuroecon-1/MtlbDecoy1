function y=nigpdf(x,mu,alpha,beta,delta)
%NIGPDF Normal-Inverse Gaussian (NIG) probability density function
%   Y = NIGPDF(X,MU,ALPHA,BETA,DELTA) returns the probability density
%   function of the Normal-Inverse Gaussian (NIG) Distribution with
%   location MU, tail heavyness ALPHA, asymmetry parameter BETA,
%   and scale parameter DELTA, at the values in X.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%        ALPHA>BETA
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 5
    error('nigpdf:TooFewInputs',...
        'Requires five input argument.');
end


[errorcode, x,mu,alpha,beta,delta] = distchck(5,x,mu,alpha,beta,delta);

if errorcode > 0
    error('nigpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(mu,'single') || isa(alpha,'single') || isa(beta,'single') || isa(delta,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

k=(alpha>beta);
if any(k)
    gamma=sqrt(alpha(k).^2-beta(k).^2);
    temp=sqrt(delta(k).^2+(x(k)-mu(k)).^2);
    term1=alpha(k).*delta(k).*besselk(1,alpha(k).*temp)./(pi*temp);
    term2=exp(delta(k).*gamma+beta(k).*(x(k)-mu(k)));
    y(k)=term1.*term2;
end

% Return NaN for out of range parameters.
y(alpha<beta)=NaN;

%Round off
y(y<0)=0;

end