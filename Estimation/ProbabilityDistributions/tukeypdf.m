function y = tukeypdf(x,lambda,mu,sigma)
%TUKEYPDF Probability density function of the Tukey Lambda Distribution
%   Y = TUKEYPDF(X,LAMBDA,MU,SIGMA) returns the probability density function
%   of the Tukey Lambda Distribution with shape parameter LAMBDA, location
%   parameter MU, and scale parameter SIGMA, at the values in X.
%
%   TUKEYPDF(X, LAMBDA) uses the default values for MU=0, SIGMA=1
%   TUKEYPDF(X, LAMBDA, MU) uses the default value SIGMA=1;
%
%   Type: Continuous, unbounded / bounded
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%


%   Mike Sheppard
%   Last Modified 31-Jul-2011




if nargin < 2
    error('tukeypdf:TooFewInputs',...
          'Requires at least two input arguments.');
end

if nargin==2
    mu=0;
    sigma=1;
end

if nargin==3
    sigma=1;
end

[errorcode, x,lambda,mu, sigma] = distchck(4,x,lambda,mu,sigma);

if errorcode > 0
    error('tukeypdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L=lambda; u=mu; s=sigma;
sz=size(L);

%Only the inverse is given
%CDF uses newton's method to find answer
%PDF needs 1/derivative for answer

%First, use CDF to find P
p = tukeycdf(x,lambda,mu,sigma);

%Now with p. calculate derivative of the inverse at those points
dinv=tukeyinvderv(p,L,u,s);

%PDF is 1/dinv;
y=1./dinv;

%Reshape
y=reshape(y,sz);

end


function dinv=tukeyinvderv(p,L,u,s)
%The parameter u is never used in the derivative
%Parameters have passed error checking
%Five cases
k1=(p>1)|(p<0);
k2=(p==0&L>0)|(p==1&L>0);
k3=(p>0&p<1&L>0)|(p>0&p<1&L<0);
k4=(L<=0&p==0)|(L<=0&p==1);
k5=~(k1|k2|k3|k4);
if any(k1)
    dinv(k1)=NaN;
elseif any(k2)
    dinv(k2)=0;
elseif any(k3)
    Lk=L(k3); pk=p(k3); sk=s(k3);
    dinv(k3)=sk.*((pk.^(Lk-1))+((1-pk).^(Lk-1)));
elseif any(k4)
    dinv(k4)=0;
else
    %k5 (L==0)
    sk=s(k5); pk=p(k5); %Note: L and u are not used
    dinv(k5)=(sk./pk)+(sk./(1-pk));
end

dinv=dinv(:);
dinv=reshape(dinv,size(p));
end