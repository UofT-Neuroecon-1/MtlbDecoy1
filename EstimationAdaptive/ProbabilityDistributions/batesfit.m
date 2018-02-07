function nhat = batesfit(x,a,b)
%BATESFIT Parameter estimate for the Bates Distribution
%   BATESFIT(X) Returns the estimate of the parameter N of the Bates
%   distribution given the data X. 
%
%   The Bates Distribution is the mean of N independent random variables
%   uniformly distributed from A to B.
%
%   Default values for A=0 and B=1.
%
%   As N is an integer, the closest integer that matches the variance will
%   be selected. In this regard it is a Method of Moments estimator.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last modified 7-May-2011


if nargin < 1
    error('batesfit:TooFewInputs',...
          'Requires at least one input argument.'); 
end

if nargin==1
  a=0;
  b=1;
end

v=var(x(:));
n = ((b-a).^2)./(12.*v);

%Non integer, check above and below for closest match to estimated variance
n = floor(n)-2:floor(n)+2; n(n==0)=[];
v2=((b-a).^2) ./ (12*n);
[val,indx] = min((v-v2).^2);
nhat=n(indx);

end