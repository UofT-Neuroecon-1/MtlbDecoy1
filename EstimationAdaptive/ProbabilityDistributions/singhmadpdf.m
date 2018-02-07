function y = singhmadpdf(x,q,a,b)
%SINGHMADPDF Singh-Maddala probability density function
%   Y = SINGHMADPDF(X,Q,A,B) returns the probability density function of
%   the Singh-Maddala Distribution with shape parameters Q and A and scale
%   parameter B, at the values in X.
%
%   This is also known as:
%        Burr XII distribution.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        X, A, B, Q > 0
%
%   Note: The Singh-Maddala distribution is a special case of the 
%   Generalized Beta Prime Distribution
%   SINGHMADPF(X,Q,A,B) = GBETAPRPDF(x,1,Q,A,B)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011


% Let gbetaprdpdf handle any errors
y = gbetaprpdf(x,1,q,a,b);



end