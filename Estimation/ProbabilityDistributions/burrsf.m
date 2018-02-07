function s = burrsf(x,c,k)
%BURRSF Burr survival function
%   S = BURRSF(X,C,K) returns the survival function of the 
%   Burr distribution with shape parameters C and K, at the values in X.  
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        C, K > 0
%
%   Note: The Burr Distribution is also known as the Burr Type XII
%   Distribution, and is a special case of the Sing-Maddala Distribution.
%   The notation for the distribution used here is F(X<x) = 1 - (1+x^C)^-K   
%
%   See also BURRPDF, BURRCDF, BURRINV, BURRSTAT, BURRFIT, BURRLIKE,
%            BURRRND, BURRHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('burrsf:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    s=1-burrcdf(x,c,k);
catch
    error('burrsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end