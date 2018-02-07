function [m,v] = singhmadstat(q,a,b)
%SINGHMADSTAT Mean and variance of the Singh-Maddala distribution
%   [M,V] = SINGHMADSTAT(Q,A,B) returns the mean and variance of the
%   Singh-Maddala distribution with shape parameters Q and A and scale
%   parameter B.

%   The Singh-Maddala distribution is a special case of the Generalized
%   Beta Prime Distribution
%   SINGHMASTAT(Q,A,B) = GBETAPRSTAT(1,Q,A,B)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011


% Let gbetaprstat handle any errors
[m,v] = gbetaprstat(1,q,a,b);


end