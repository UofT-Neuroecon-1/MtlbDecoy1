function [m,v] = baldnichstat(F,p)
%BALDNICHSTAT Mean and variance for the Balding-Nichols distribution
%   [M,V] = BALDNICHSTAT(F,P) returns the mean and variance of the
%   Balding-Nichols Distribution with parameters F and p.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
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
%   See also BALDNICHPDF, BALDNICHCDF, BALDNICHINV, BALDNICHFIT,
%            BALDNICHLIKE, BALDNICHRND, BALDNICHSF, BALDNICHHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin < 2
    error('baldnichstat:TooFewInputs',...
        'Requires at least two input arguments.');
end



try
    %Parametrization of Beta
    [m,v] = betastat((1-F).*p./F,(1-F).*(1-p)./F);
catch err
    error('baldnichstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<F & F<1) & (0<p & p<1);
m(~okparam)=NaN;
v(~okparam)=NaN;


end