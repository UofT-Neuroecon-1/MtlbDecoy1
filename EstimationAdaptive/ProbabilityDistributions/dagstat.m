function [m,v] = dagstat(p,a,b)
%DAGSTAT Mean and variance for the Dagum distribution
%   [M,V] = DAGSTAT(P,A,B) returns the mean and variance of the Dagum
%   distribution with shape parameters P and A and scale parameter B.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        P, A , B > 0
%
%   Note: The Dagum Distribution is also known as the Burr III
%   Distribution, and is also a special case of the Generalized 
%   Beta Prime Distribution.
%
%   See also DAGPDF, DAGCDF, DAGINV, DAGFIT, DAGLIKE, DAGRND, DAGSF, DAGHAZ
%

%   Mike Sheppard
%   Last Modified 7-May-2011

if nargin < 3
    error('dagstat:TooFewInputs',...
          'Requires at least three input argument.'); 
end

%Special case of the Generalized Beta Prime Distribution
[m,v]=gbetaprstat(p,1,a,b);

end
