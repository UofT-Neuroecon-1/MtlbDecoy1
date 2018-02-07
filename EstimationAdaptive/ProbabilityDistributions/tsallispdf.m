function y = tsallispdf(x,q)
%TSALLISPDF Probability density function for the Tsallis Distribution
%   Y = TSALLISPDF(X,Q) returns the probability density function of the
%   Tsallis distribution with Q degrees of freedom, at the values in X.
%
%   The Tsallis distribution is a reparameterization of Student's T
%   distribution with V degrees of freedom in the T distribution given by
%   T = (3-Q) / (Q-1)
%
%   Type: Continuous, unbounded
%   Restrictions:
%        1 < Q < 2
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Dec-2011


if nargin < 2,
    error('tsallispdf:TooFewInputs');
end

% Return NaN for out of range parameters.
q(q<1 | q>2)=NaN;

try
    v=(3-q)./(q-1);
    %Use Student's T distribution to solve for answer and error check
    y = tpdf(x,v);
catch
    error('tsallispdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end

