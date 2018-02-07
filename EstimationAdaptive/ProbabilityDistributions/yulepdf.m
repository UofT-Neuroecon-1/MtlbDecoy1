function y = yulepdf(x,a)
%YULEPDF Yule–Simon probability density function
%   Y = YULEPDF(X,A) returns the Yule-Simon density function with
%   shape parameter A, at the values in X.
%
%   Note: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 30-May-2011



if nargin < 2
    error('yulepdf:TooFewInputs',...
        'Requires at least two input arguments.');
end


try
    %Let Beta catch errors on x and a, yielding NaNs as outputs
    y = a.*beta(1+x,a+1);
catch
    error('yulepdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end