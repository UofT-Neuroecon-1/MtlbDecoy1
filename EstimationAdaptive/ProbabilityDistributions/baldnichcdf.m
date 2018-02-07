function y = baldnichcdf(x,F,p)
%BALDNICHCDF Balding–Nichols cumulative distribution function
%   Y = BALDNICHCDF(x,F,P) returns the cumulative distribution function
%   of the Balding–Nichols Distribution with parameters F and p, at the
%   values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The distribution assumes background allele frequency P the allele
%   frequencies, in sub-populations separated by Wright's F_ST F.
%
%   Distribution: Continuous, bounded, (0,1)
%   Restrictions:
%        0 < F, P < 1
%
%   Note: The Balding-Nichols distribution is a reparametrization of the
%   Beta distribution
%
%   See also BALDNICHPDF, BALDNICHINV, BALDNICHSTAT, BALDNICHFIT, 
%            BALDNICHLIKE, BALDNICHRND, BALDNICHSF, BALDNICHHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin ~= 3
    error('baldnichcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x F p] = distchck(3,x,F,p);

if errorcode > 0
    error('baldnichpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<F & F<1) & (0<p & p<1);

%Parametrization of Beta
y = betacdf(x,(1-F).*p./F,(1-F).*(1-p)./F);

y(okparam & x<0)=0;
y(okparam & x>1)=1;
y(~okparam)=NaN;

%Catch round off
y(y<0)=0; y(y>1)=1;

end