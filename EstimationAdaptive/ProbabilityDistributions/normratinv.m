function x = normratinv(p,s1,s2)
%NORMRATINV Inverse of the Normal Ratio cumulative distribution function
%   X = NORMRATINV(P,S1,S2) returns the inverse of the Normal Ratio
%   cumulative distribution function of the ratio of two Normal
%   Distributions with standard deviations S1 and S2.
%
%   This is equivalent to the Cauchy Distribution
%   N(0,S1)/N(0,S2)=Cauchy(S1/S2)
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      s1 , s2 >0
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 3
    error('normratinv:TooFewInputs',...
          'Requires at three input argument.');
end

[errorcode, p,s1,s2] = distchck(3,p,s1,s2);

if errorcode > 0
    error('normratinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

scale=s1./s2;
x = cauchyinv(p,0,scale);

end