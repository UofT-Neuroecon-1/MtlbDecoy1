function parmhat = asinfit(x)
%ASINFIT Parameter estimates for the Arcsine Distribution
%   PARMHAT = ASINFIT(X) returns maximum likelihood estimates of the
%   parameters of a Arcsine Distribution fit to the data in X. PARMHAT(1)
%   and PARMHAT(2) are estimates of A and B, respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINPDF, ASINCDF, ASININV, ASINSTAT, 
%            ASINLIKE, ASINRND, ASINSF, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 25-May-2012


if nargin < 1
    error('asinfit:TooFewInputs','Requires at least one input argument.'); 
end

x=x(:);
parmhat = [min(x) max(x)];

end
