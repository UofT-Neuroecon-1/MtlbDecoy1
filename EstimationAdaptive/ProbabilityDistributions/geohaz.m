function h = geohaz(x,p)
%GEOHAZ Geometric hazard function.
%   H = GEOHAZ(X,P) returns the hazard function of the geometric 
%   distribution with probability parameter P, evaluated at the 
%   values in X. 
%
%   The size of H is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   See also GEOPDF, GEOCDF, GEOINV, GEOSTAT, GEOFIT, GEOLIKE, 
%            GEORND, GEOSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 2
    error(message('geohaz:TooFewInputs'));
end

try
    yt = geopdf(x,p);
    st = geosf(x,p);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('geohaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
