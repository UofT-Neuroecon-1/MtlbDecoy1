function p = loglappdf(x,u,b)
%LOGLAPCDF Log-Laplace probability density function
%   y = LOGLAPPDF(X,U,B) returns the Log-Laplace probability density
%   function with distribution parameters U and B. U and B are the mean and
%   scale parameter, respectively, of the associated Laplace distribution.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     X>0
%     B>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('loglappdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,u,b] = distchck(3,x,u,b);

if errorcode > 0
    error('loglappdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Transformation of the Laplace Distribution
y = lappdf(log(x),u,b)./x;

%Valid only for x>0
y(x<0)=0;

end