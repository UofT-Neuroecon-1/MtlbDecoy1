function [m,v]=bnkgbstat(a,b)
%BNKGBSTAT Mean and variance for the Benktander-Gibrat distribution
%   [M,V] = BNKGBSTAT(A,B) returns the mean and variance of the
%   Benktander-Gibrat distribution with parameters A and B.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%        A(A+1) >= 2B
%        A , B > 0
%
%   Note: The Benktander-Gibrat Distribution is also known as the
%   Benktander Distribution of Type I.
%
%   See also BNKGBPDF, BNKGBCDF, BNKGBINV, BNKGBFIT, BNKGBLIKE, BNKGBRND,
%            BNKGBSF, BNKGBHAZ
%

%   Mike Sheppard
%   Last Modified 13-Mar-2011


if nargin < 2
    error('bnkgbstat:TooFewInputs',...
        'Requires at least two input arguments.');
end



try
    %Expand size if necessary
    a=a+zeros(size(b));
    b=b+zeros(size(a));
    
    % Return NaN if any arguments are outside of their respective limits.
    okparam = (0<a & a<Inf) & (0<b & b<Inf) & (a.*(a+1) >= 2*b);

    %Mean
    m=1+(1./a);
    
    %variance compute in terms
    term1=(-1+a)./(2.*sqrt(b));
    term2=-sqrt(b);
    term3=a.*exp(term1.^2);
    term4=erfc(term1);
    term5=(a.^2).*sqrt(b);
    v=(term2+(term3.*sqrt(pi).*term4))./term5;
    
    m(~okparam)=NaN;
    v(~okparam)=NaN;
    
catch
    error('bnkgbstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end



end