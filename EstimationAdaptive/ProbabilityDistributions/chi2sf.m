function s = chi2sf(x,v)
%CHI2SF Chi-square survival function.
%   S = CHI2SF(X,V) returns the survival function of the chi-square 
%   distribution with V degrees of freedom at the values in X.
%   The chi-square density function with V degrees of freedom,
%   is the same as a gamma density function with parameters V/2 and 2.
%
%   The size of S is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.      
%
%   See also CHI2PDF, CHI2CDF, CHI2INV, CHI2STAT, 
%            CHI2FIT, CHI2LIKE, CHI2RND, CHI2HAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 2
    error('chi2sf:TooFewInputs',...
        'Requires two input argument.');
end


try
    s = 1 - chi2cdf(x,v);
catch
    error('chi2sf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
