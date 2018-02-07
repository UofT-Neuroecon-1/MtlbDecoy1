function x = daginv(y,p,a,b)
%DAGINV Inverse of the Dagum cumulative distribution function
%   X = DAGINV(Y,P,A,B) returns the inverse cumulative distribution
%   function of the Dagum Distribution with shape parameters P and A
%   and scale parameter B, at the values in Y.
%
%   The size of X is the common sizes of the inputs. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%      P, A , B > 0
%
%   NOTE: The Dagum Distribution is also known as the Burr III
%   Distribution, and is also a special case of the Generalized
%   Beta Prime Distribution.
%
%   See also DAGPDF, DAGCDF, DAGSTAT, DAGFIT, 
%            DAGLIKE, DAGRND, DAGSF, DAGHAZ
%

%   Mike Sheppard
%   Last Modified 7-May-2011

if nargin ~= 4
    error('daginv:TooFewInputs',...
        'Requires four input arguments.');
end

try
    %Special case of the Generalized Beta Prime Distribution
    x=gbetaprinv(y,p,1,a,b);
catch err
    error('daginv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end