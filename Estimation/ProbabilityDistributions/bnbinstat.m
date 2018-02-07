function [m,v] = bnbinstat(n,a,b)
%BNBINSTAT Mean and variance for the Beta Negative Binomial Distribution
%   [M,V] = BNBINSTAT(N,A,B) returns the mean and variance of the Beta
%   Negative Binomial Distribution with parameters A and B, 
%   with N successes.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Discrete, semi-bounded, {0,...,Inf}
%   Restrictions:
%        A , B > 0
%            N >= 0 (integer)
%
%   See also BNBINPDF, BNBINCDF, BNBININV, BNBINFIT, BNBINLIKE, BNBINRND,
%            BNBINSF, BNBINHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011


if nargin < 3
    error('bbinostat:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode n a b] = distchck(3,n,a,b);

if errorcode > 0
    error('bbinostat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=0 & n==round(n));

m=n.*b./(a-1);
v=(n.*(a+n-1).*b.*(a+b-1))./((a-2).*((a-1).^2));

m(~okparam)=NaN; v(~okparam)=NaN;
m(okparam & a<=1)=Inf;
v(okparam & a<=2)=Inf;


end