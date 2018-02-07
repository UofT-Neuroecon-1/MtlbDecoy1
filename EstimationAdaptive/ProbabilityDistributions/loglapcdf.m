function p = loglapcdf(x,u,b)
%LOGLAPCDF Log-Laplace cumulative distribution function
%   P = LOGLAPCDF(X,U,B) returns the Log-Laplace cumulative distribution
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
   error('loglapcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,u,b] = distchck(3,x,u,b);

if errorcode > 0
    error('loglapcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Transformation of the Laplace Distribution
p = lapcdf(log(x),u,b);

%Valid only for x>0
p(x<0)=0;

end