function x = loglapinv(p,u,b)
%LOGLAPINV Inverse of the Log-Laplace cumulative distribution function
%   X = LOGLAPINV(P,U,B) returns the Log-Laplace cumulative distribution
%   function with distribution parameters U and B. U and B are the mean and
%   scale parameter, respectively, of the associated Laplace distribution.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     0<=P<=1
%     B>0
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('loglapinv:TooFewInputs','Requires three input arguments.');
end

[errorcode, p,u,b] = distchck(3,p,u,b);

if errorcode > 0
    error('loglapinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Transformation of the Laplace Distribution
x = exp(lapinv(p,u,b));

end