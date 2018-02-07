function x = baldnichinv(y,F,p)
%BALDNICHINV Inverse of the Balding–Nichols cumulative distribution function
%   X = BALDNICHINV(Y,F,P) returns the inverse cumulative distribution
%   function of the Balding-Nichols Distribution with parameters F and p,
%   at the values in y.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The distribution assumes background allele frequency P the allele
%   frequencies, in sub-populations separated by Wright's F_ST F.
%
%   Distribution: Continuous, bounded, (0,1)
%   Restrictions:
%      0 < F, P < 1
%
%   NOTE: The Balding-Nichols distribution is a reparametrization of the
%   Beta distribution
%
%   See also BALDNICHPDF, BALDNICHCDF, BALDNICHSTAT, BALDNICHFIT,
%            BALDNICHLIKE, BALDNICHRND, BALDNICHSF, BALDNICHHAZ
%

%   Mike Sheppard
%   Last Modified 16-May-2012


if nargin ~= 3
    error('baldnichinv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode y F p] = distchck(3,y,F,p);

if errorcode > 0
    error('baldnichinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize X to zero.
if isa(y,'single') || isa(F,'single') || isa(p,'single')
    x = zeros(size(y),'single');
else
    x = zeros(size(y));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<F & F<1) & (0<p & p<1);
okvar = (0 < y & y < 1);
ok = (okparam & okvar);
x(~ok)=NaN;
%Edge cases
x(okparam & y==0)=0;
x(okparam & y==1)=1;


if any(ok)
    y=y(ok); F=F(ok); p=p(ok);
    %Parametrization of Beta
    x(ok) = betainv(y,(1-F).*p./F,(1-F).*(1-p)./F);
end



end
