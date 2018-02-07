function s = geosf(x,p)
%GEOSF Geometric survival function.
%   S = GEOSF(X,P) returns the survival function of the geometric 
%   distribution with probability parameter P, evaluated at the 
%   values in X. 
%
%   The size of S is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   See also GEOPDF, GEOCDF, GEOINV, GEOSTAT, GEOFIT, GEOLIKE, 
%            GEORND, GEOHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 2
    error(message('geosf:TooFewInputs'));
end

try
    s = 1 - geocdf(x,p);
catch
    error('geosf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

