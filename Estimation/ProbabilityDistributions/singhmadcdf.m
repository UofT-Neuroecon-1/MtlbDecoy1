function y = singhmadcdf(x,q,a,b)
%SINGHMADCDF Singh-Maddala cumulative distribution function
%   Y = SINGHMADCDF(X,Q,A,B) returns the Singh-Maddala cumulative distribution
%   function with shape parameters Q and A and scale parameter B at the
%   values in X.
%
%   The Singh-Maddala distribution is a special case of the Generalized
%   Beta Prime Distribution
%   SINGHMACPF(X,Q,A,B) = GBETAPRCDF(x,1,Q,A,B)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011


% Let gbetaprcdf handle any errors
y = gbetaprcdf(x,1,q,a,b);


end