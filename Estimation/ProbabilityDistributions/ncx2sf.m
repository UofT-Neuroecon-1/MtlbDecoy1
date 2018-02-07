function s = ncx2sf(x,v,delta)
%NCX2SF Non-central chi-square survival function
%   S = NCX2SF(X,V,DELTA) Returns the survival function of the non-central
%   chi-square distribution with V degrees of freedom and non-centrality
%   parameter, DELTA, at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Some texts refer to this distribution as the generalized Rayleigh,
%   Rayleigh-Rice, or Rice distribution.
%
%   See also NCX2PDF, NCX2CDF, NCX2INV, NCX2STAT, NCX2FIT, NCX2LIKE,
%            NCX2RND, NCX2HAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011



if nargin < 3
    error(message('ncx2sf:TooFewInputs'));
end


try
    s = 1 - ncx2cdf(x,v,delta);
catch
    error('ncx2sf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
