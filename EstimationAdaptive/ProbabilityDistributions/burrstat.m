function [m,v]=burrstat(c,k)
%BURRSTAT Mean and variance for the Burr distribution
%   [M,V] = BURRSTAT(C,K) returns the mean and variance of the Burr
%   distribution with shape parameters C and K.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        C, K > 0
%
%   Note: The Burr Distribution is also known as the Burr Type XII
%   Distribution, and is a special case of the Sing-Maddala Distribution.
%   The notation for the distribution used here is F(X<x) = 1 - (1+x^C)^-K
%
%   See also BURRPDF, BURRCDF, BURRINV, BURRFIT, BURRLIKE, BURRRND, 
%            BURRSF, BURRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011



if nargin < 2
    error('burrstat:TooFewInputs',...
        'Requires at least two input arguments.');
end


try
    %Moments uses the gammaln for accuracy
    m1=exp(gammaln((1+c)./c)+gammaln(k-(1./c))-gammaln(k));
    m2=exp(gammaln((2+c)./c)+gammaln(k-(2./c))-gammaln(k));
    
    % Return NaN if any arguments are outside of their respective limits.
    okparam = (0<c & c<Inf) & (0<k & k<Inf);
    
    m=m1;
    v=m2-(m1.^2);
    
    m(~okparam)=NaN; v(~okparam)=NaN;
    
catch
    error('burrstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any of A,B or S are not positive.
m(c<0 | k<0 ) = NaN;
v(c<0 | k<0 ) = NaN;

end