function [m,v] = betaprstat(a,b)
%BETAPRSTAT Mean and variance for the Beta Prime Distribution.
%   [M,V] = BETAPRSTAT(A,B) returns the mean and variance of the
%   Beta Prime distribution with shape parameters A and B
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRPDF, BETAPRCDF, BETAPRINV, BETAPRFIT,
%            BETAPRLIKE, BETAPRRND, BETAPRSF, BETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin < 2
    error('betaprstat:TooFewInputs',...
        'Requires at least two input arguments.');
end


try
    %Expand size if necessary
    a=a+zeros(size(b));
    b=b+zeros(size(a));
    % Return NaN if any arguments are outside of their respective limits.
    okparam = (0<a & a<Inf) & (0<b & b<Inf);
    
    m=a./(b-1);
    v=a.*(b+a-1)./((b-2).*((b-1).^2));
    
    m(~okparam)=NaN; m(okparam & b<=1)=Inf;
    v(~okparam)=NaN; v(okparam & b<=2)=NaN;
catch
    error('betaprstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end