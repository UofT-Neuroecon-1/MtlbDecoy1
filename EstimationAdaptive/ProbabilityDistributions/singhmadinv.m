function x = singhmadinv(y,q,a,b)
%SINGHMADINV Inverse of the Singh-Maddala CDF
%   X = SINGHMADINV(Y,Q,A,B) returns the inverse cumulative distribution 
%   function of the Singh-Maddala distribution with shape parameters Q 
%   and A and scale parameter B at the values in Y.
%
%   The Singh-Maddala distribution is a special case of the Generalized
%   Beta Prime Distribution
%   SINGHMAINV(X,Q,A,B) = GBETAPRINV(x,1,Q,A,B)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011


% Let gbetaprinv handle any errors
x = gbetaprinv(y,1,q,a,b);


end