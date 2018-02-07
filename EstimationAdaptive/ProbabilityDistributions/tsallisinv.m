function x = tsallisinv(p,q)
%TSALLISINV Inverse for Tsallis distribution
%   X = TSALLISINV(P,Q) returns the inverse  of the Tsallis distribution with Q
%   degrees of freedom, at the values in P.
%
%   The Tsallis distribution is just a reparameterization of Student's T
%   distribution
%   
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.   
%

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 2, 
    error('tsalliscdf:TooFewInputs'); 
end

v=(3-q)./(q-1);

%Use Student's T distribution to solve for answer and error check
x = tinv(p,v);

end

