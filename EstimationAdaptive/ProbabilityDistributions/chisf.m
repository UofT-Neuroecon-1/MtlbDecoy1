function s = chisf(x,v)
%CHISF Chi survival function.
%   S = CHISF(X,V) returns the survival function of the Chi
%   Distribution with V degrees of freedom, at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        V > 0
%
%   See also CHIPDF, CHICDF, CHIINV, CHISTAT, CHIFIT, CHILIKE, 
%            CHIRND, CHIHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if   nargin < 2,
    error('chisf:TooFewInputs','Requires two input arguments.');
end


try
    s=1-chicdf(x,v);
catch
    error('chisf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end
