function [m,v] = normratstat(s1,s2)
%NORMRATSTAT Mean and variance of Normal Ratio Distribution
%   X = NORMRATSTAT(S1,S2) returns the mean and variance of the Normal
%   Ratio Distribution, of the ratio of two Normal
%   Distributions with standard deviations S1 and S2.
%
%   This is equivalent to the Cauchy Distribution
%   N(0,S1)/N(0,S2)=Cauchy(S1/S2)
%
%   As it is identical to the Cauchy Distribution the mean and variance do
%   not exist
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 2
    error('normsumstat:TooFewInputs',...
        'Requires at least two input argument.');
end


try
    scale=s1./s2;
    [m,v] = cauchystat(0,scale);
catch
    error('normratstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end