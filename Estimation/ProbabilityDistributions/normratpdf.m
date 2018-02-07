function y = normratpdf(x,s1,s2)
%NORMRATPDF Normal Ratio Distribution probability density function
%   Y = NORMRATPDF(X,S1,S2) returns the probability density function of
%   the Normal Ratio Distribution of the ratio of two Normal Distributions
%   with standard deviations S1 and S2.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      s1 , s2 >0
%
%   Note: This is equivalent to the Cauchy Distribution with
%   parameter S1/S2
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 3
    error('normratpdf:TooFewInputs',...
          'Requires at three input argument.');
end

[errorcode, x,s1,s2] = distchck(3,x,s1,s2);

if errorcode > 0
    error('normratpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN for out of range parameters.
s1(s1<0)=NaN; s2(s2<0)=NaN;

scale=s1./s2;
y = cauchypdf(x,0,scale);

%Round off
y(y<0)=0;

end