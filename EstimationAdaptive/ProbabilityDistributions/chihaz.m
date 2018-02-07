function h = chihaz(x,v)
%CHIHAZ Chi hazard function.
%   H = CHIHAZ(X,V) returns the hazard function of the Chi
%   Distribution with V degrees of freedom, at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        V > 0
%
%   See also CHIPDF, CHICDF, CHIINV, CHISTAT, CHIFIT, CHILIKE, 
%            CHIRND, CHISF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if  nargin < 2,
    error('chihaz:TooFewInputs','Requires two input arguments.');
end


try
    h = chipdf(x,v) ./ chisf(x,v);
catch
    error('chihaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end
