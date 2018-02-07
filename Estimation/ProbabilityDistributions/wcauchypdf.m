function y = wcauchypdf(theta,mu,lambda)
%WCAUCHYPDF Wrapped Cauchy probability density function
%   Y = WCAUCHYPDF(theta,mu,lambda) returns the Wrapped Cauchy probability density
%   function with location parameter MU and scale parameter LAMBDA.
%
%   Default values are MU=0 and LAMBDA=1.
%
%   The size of Y is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%

%   Mike Sheppard
%   Last Modified 6-Apr-2011

if nargin<1
    error('wcauchypdf:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2
    mu=0;
end
if nargin < 3
    lambda=1;
end


[errorcode theta mu lambda] = distchck(3,theta,mu,lambda);

if errorcode > 0
    error('wcauchypdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(theta,'single') || isa(mu,'single') || isa(lambda,'single')
    y=zeros(size(theta),'single');
else
    y=zeros(size(theta));
end


% Return NaN if A is not real or B not real positive
y( ~isreal(mu) | ~isreal(mu) | lambda<0 ) = NaN;

k1=find( isreal(mu) & isreal(mu) & isreal(lambda) & lambda>0 );
if any(k1),
    thetak=theta(k1);
    muk=mu(k1);
    lambdak=lambda(k1);
    y(k1)=(sinh(lambdak)./(cosh(lambdak)-cos(theta-mu)))./(2*pi);
end

end