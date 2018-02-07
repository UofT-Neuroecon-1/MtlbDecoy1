function x = loginv(p,mu,sigma)
%LOGINV Inverse Logistic cumulative distribution function
%   X = LOGINV(P,MU,SIGMA) returns the inverse logistic cumulative function
%   with mean U and scale parameter SIGMA.
%
%   NOTE: loginv(p,0,1) is equivalent to logit(p)
%
%   Type: Continuous, unbounded
%   Restrictions:
%     0<=P<=1
%     SIGMA>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('loginv:TooFewInputs','Requires three input arguments.');
end

[errorcode, p,mu,sigma] = distchck(3,p,mu,sigma);

if errorcode > 0
    error('loginv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

k=(sigma>0);
if any(k)
    x(k) = mu(k) + sigma(k).*logit(p(k));  %LOGIT function below
end

%Else not valid
x(~k)=NaN;


end



function lp=logit(p)
lp=zeros(size(p));
%Compute it accurately using log1p function
%Logit(p) == log(p) - log(1-p);

k1=(p<.1); %(1-p) is close to 1
if any(k1)
    lp(k1)=log(p(k1))-log1p(-p(k1)); %log(p)-log(1+(-p))
end

k2=(p>.9);  %p is close to 1
if any(k2)
    lp(k2)=log1p(p(k2)-1)-log(1-p(k2)); %log(1+(p-1))-log(1-p)
end

k3=(~k1)&(~k2);
if any(k3)
    lp(k3)=log(p(k3))-log(1-p(k3));
end
end
