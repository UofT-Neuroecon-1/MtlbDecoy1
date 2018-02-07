function y = tsalliscdf(x,q)
%TSALLISCDF cumulative distribution function (cdf) for Tsallis distribution
%   Y = TSALLISCDF(X,Q) returns the pdf of the Tsallis distribution with Q
%   degrees of freedom, at the values in X.
%
%   The Tsallis distribution is just a reparameterization of Student's T
%   distribution
%   
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 2, 
    error('tsalliscdf:TooFewInputs'); 
end

v=(3-q)./(q-1);

%Use Student's T distribution to solve for answer and error check
y = tcdf(x,v);

end

