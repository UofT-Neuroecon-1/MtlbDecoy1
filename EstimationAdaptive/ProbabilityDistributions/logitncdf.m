function p = logitncdf(x,mu,sigma)
%LOGITNCDF Logit-normal cumulative distribution function
%   P = LOGITNCDF(X,MU,SIGMA) returns the cumulative distribution 
%   function of the Logit-normal distribution with parameters MU and SIGMA. 
%   MU and SIGMA are the mean and  standard deviation, respectively, 
%   of the associated normal distribution.  
%
%   Type: Continuous, bounded
%   Restrictions:
%     0<=X<=1
%     SIGMA>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('logitncdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('logitncdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

k=(x>=0 & x<=1 & sigma>0);
if any(k)
%Transformation of the Normal Distribution
%Compute logit from Logistic Inverse
logitx=loginv(x(k),0,1);
p(k) = lapcdf(logitx,mu(k),sigma(k));
end

%Else not valid
p(~k)=NaN;

end
