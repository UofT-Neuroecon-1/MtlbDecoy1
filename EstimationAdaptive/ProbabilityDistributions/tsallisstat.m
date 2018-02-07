function [m,v]=tsallisstat(q)
%TSALLISSTAT Mean and variance for the Tsallis distribution
%   [M,V] = TSALLISSTAT(Q) returns the mean and variance of Tsallis
%   Distribution with Q degrees of freedom
%
%   The Tsallis distribution is just a reparameterization of Student's T
%   distribution
%   
%   The size of out arguments is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.   
%

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 1, 
    error('tsallisstat:TooFewInputs'); 
end

v=(3-q)./(q-1);

%Use Student's T distribution to solve for answer and error check
[m,v]=tstat(v);

end
