function [m,v] = rbeckstat(s1,s2)
%RBECKSTAT Mean and variance for the Reduced Beckmann distribution
%   [M,V] = RBECKSTAT(S1,S2) returns the mean and variance of the
%   Reduced Beckmann distribution with standard deviations S1 and S2.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   The full Beckmann Distribution
%   If {x,y} follows the Bivariate Normal Distribution with means {u1,u2},
%   standard deviations {s1,s2} and correlation rho, then sqrt[x^2+y^2]
%   follows the Beckmann Distribution [u1,u2,s1,s2,rho]
%
%   The reduced Beckmann distribution assumes u1=u2=rho=0 so the bivariate
%   normal distributions are centered at (0,0) with covariance matrix
%   [s1^2 0; 0 s2^2]
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        S1, S2 > 0
%
%   See also RBECKPDF, RBECKCDF, RBECKINV, RBECKFIT, RBECKLIKE, RBECKRND, 
%            RBECKSF, RBECKHAZ
%

%   Mike Sheppard
%   Last Modified 19-Dec-2011


if nargin < 2
    error('rbeckstat:TooFewInputs',...
        'Requires at least two input arguments.');
end

try
    %Expand size if necessary
    s1=s1+zeros(size(s2));
    s2=s2+zeros(size(s1));
    
    %Make NaN for invalid inputs
    s1(s1<0) = NaN;
    s2(s2<0) = NaN;
 
    b=1-(min(s1,s2)./max(s1,s2)).^2;
    [K,E] = ellipke(b);
    m=max(s1,s2).*(sqrt(2/pi)).*E;
    v=(s1.^2+s2.^2)-(m.^2);
    
catch
    error('rbeckstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end



end