function y = baldnichpdf(x,F,p)
%BALDNICHPDF Balding–Nichols probability density function
%   Y = BALDNICHPDF(x,F,P) returns the probability density function of the
%   Balding–Nichols distribution with parameters F and p at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The distribution assumes background allele frequency P the allele
%   frequencies, in sub-populations separated by Wright's F_ST F.
%
%   Distribution: Continuous, bounded, (0,1)
%   Restrictions:
%        0 < F , P < 1
%
%   Note: The Balding-Nichols distribution is a reparametrization of the
%   Beta distribution
%
%   See also BALDNICHCDF, BALDNICHINV, BALDNICHSTAT, BALDNICHFIT, 
%            BALDNICHLIKE, BALDNICHRND, BALDNICHSF, BALDNICHHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin < 3
    error('baldnichpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x F p] = distchck(3,x,F,p);

if errorcode > 0
    error('baldnichpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<F & F<1) & (0<p & p<1);
okvar = (0 < x & x < 1);

%Parametrization of Beta, let Beta catch most errors
y = betapdf(x,(1-F).*p./F,(1-F).*(1-p)./F);

y(okparam & ~okvar)=0;
y(~okparam)=NaN;

%Catch round off
y(y<0)=0;

end