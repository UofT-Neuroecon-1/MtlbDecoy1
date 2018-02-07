function [m,v]=bnkwblstat(a,b)
%BNKWBLSTAT Mean and variance for the Benktander-Weibull distribution
%   [M,V] = BNKWBLSTAT(A,B) returns the mean and variance of the
%   Benktander-Weibull distribution with parameters A and B.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%        A > 0
%        0 < B <= 1
%
%   Note: The Benktander-Weibull Distribution is also known as the
%   Benktander Distribution of Type II
%
%   See also BNKWBLPDF, BNKWBLCDF, BNKWBLINV, BNKWBLFIT,
%            BNKWBLLIKE, BNKWBLRND, BNKWBLSF, BNKWBLHAZ
%

%   Mike Sheppard
%   Last Modified: 16-Dec-2011


if nargin < 2
    error('bnkwblstat:TooFewInputs',...
        'Requires at least two input arguments.');
end



try
    %Expand size if necessary
    a=a+zeros(size(b));
    b=b+zeros(size(a));
    
    % Return NaN if any arguments are outside of their respective limits.
    okparam = (0<a & a<Inf) & (0<b & b<=1);
    a(~okparam)=NaN; b(~okparam)=NaN;

    m=1+(1./a);
    %Variance in terms of transformed exponential integral function En(z)
    %Generalized Exponential Integral, in form of incomplete gamma
    n=1-(1./b); x=a./b;
    gexp=(x.^(n-1)).*gamma(1-n).*gammainc(x,1-n,'upper'); %En(z)
    v=(-1+(2.*a.*exp(a./b).*gexp./b))./(a.^2);
catch
    error('bnkwblstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end