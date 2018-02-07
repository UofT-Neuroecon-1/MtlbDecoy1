function s = invgamsf(x,a,b)
%INVGAMSF Inverse-Gamma survival function
%   S = INVGAMSF(X,A,B) returns the survival function of the Inverse-Gamma
%   Distribution with shape parameter A and scale parameter B, evaluated
%   at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%     A, B > 0
%
%   See also INVGAMPDF, INVGAMCDF, INVGAMINV, INVGAMSTAT, 
%            INVGAMFIT, INVGAMLIKE, INVGAMRND, INVGAMHAZ
%

%   Mike Sheppard
%   Last Modified: 13-May-2012


if nargin < 3
   error('invgamsf:TooFewInputs','Requires three input arguments.');
end


try
    s = 1 - invgamcdf(x,a,b);
catch
    error('invgamsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end