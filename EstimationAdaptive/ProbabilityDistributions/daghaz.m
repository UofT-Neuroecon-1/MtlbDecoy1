function h = daghaz(x,p,a,b)
%DAGHAZ Dagum hazard function
%   H = DAGHAZ(X,P,A,B) returns the hazard function of the Dagum 
%   Distribution with shape parameters P and A and scale
%   parameter B, at the values in X.
%
%   The size of S is the common sizes of the inputs. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        P, A , B > 0
%
%   Note: The Dagum Distribution is also known as the Burr III
%   Distribution, and is also a special case of the Generalized
%   Beta Prime Distribution.
%
%   See also DAGPDF, DAGCDF, DAGINV, DAGSTAT, DAGFIT, DAGLIKE, 
%            DAGRND, DAGSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 4
    error('daghaz:TooFewInputs',...
        'Requires at least four input argument.');
end

try
    h = dagpdf(x,p,a,b) ./ dagsf(x,p,a,b);
catch
    error('daghaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end
