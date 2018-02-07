function y = logitnpdf(x,mu,sigma)
%LOGITNPDF Logit-normal probability density function
%   Y = LOGITNPDF(X,MU,SIGMA) returns the probability density function of
%   the Logit-Normal Distribution with parameters MU and SIGMA. 
%
%   MU and SIGMA are the mean and standard deviation, respectively, 
%   of the associated Normal Distribution.  
%
%   Type: Continuous, bounded
%   Restrictions:
%     0 <= X <= 1
%     SIGMA > 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('logitnpdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('logitnpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

k=(x>=0 & x<=1 & sigma>0);
if any(k)
%Transformation of the Normal Distribution
logitx=loginv(x(k),0,1); %Compute logit from Logistic Inverse
dlogit=1./(x(k).*(1-x(k)));
y(k) = lappdf(logitx,mu(k),sigma(k)).*dlogit;
end

%Else zero or not valid
y(sigma<0)=NaN;
y(sigma>0 & (x<0 | x>1))=0;

end
