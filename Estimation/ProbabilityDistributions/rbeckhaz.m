function h = rbeckhaz(x,s1,s2)
%RBECKHAZ Reduced Beckmann hazard function
%   H = RBECKHAZ(X,S1,S2) returns the Reduced Beckmann hazard function
%   with standard deviations S1 and S2, at the values of X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
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
%   Type: Continuous, semi-bounded, [0,Inf) 
%   Restrictions:
%        S1, S2 > 0
%
%   See also RBECKPDF, RBECKCDF, RBECKINV, RBECKSTAT, RBECKFIT, RBECKLIKE,
%            RBECKRND, RBECKSF
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 3
    error('rbeckhaz:TooFewInputs',...
        'Requires at least three input arguments.');
end


try
    h = rbeckpdf(x,s1,s2) ./ rbecksf(x,s1,s2);
catch
    error('rbeckhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end