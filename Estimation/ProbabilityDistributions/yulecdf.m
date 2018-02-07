function p = yulecdf(x,a)
%YULECDF Yule–Simon cumulative distribution function
%   P = YULECDF(X,A) returns the Yule-Simon cumulative distribution
%   function with shape parameter A, at the values in X.
%
%   Note: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 2
    error('yulecdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end


%Let Beta catch errors on x and a, yielding NaNs as outputs
try
p = 1-a.*beta(a,x+2);
catch
        error('yulecdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


end