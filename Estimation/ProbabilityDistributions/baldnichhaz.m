function h = baldnichhaz(x,F,p)
%BALDNICHHAZ Balding–Nichols hazard function
%   H = BALDNICHHAZ(x,F,P) returns the hazard function of the
%   Balding–Nichols Distribution with parameters F and p, at the
%   values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The distribution assumes background allele frequency P the allele
%   frequencies, in sub-populations separated by Wright's F_ST F.
%
%   Type: Continuous, bounded, (0,1)
%   Restrictions:
%        0 < F, P < 1
%
%   Note: The Balding-Nichols distribution is a reparametrization of the
%   Beta distribution
%
%   See also BALDNICHPDF, BALDNICHCDF, BALDNICHINV, BALDNICHSTAT, 
%            BALDNICHFIT, BALDNICHLIKE, BALDNICHRND, BALDNICHSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 3
    error('baldnichhaz:TooFewInputs',...
        'Requires at least three input arguments.');
end


try
    h = baldnichpdf(x,F,p) ./ baldnichsf(x,F,p);
catch
    error('baldnichhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
