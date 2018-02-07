function r = beckrnd(u1,u2,s1,s2,rho,varargin)
%BECKRND Random arrays from the Beckmann distribution
%   R = BECKRND(U1,U2,S1,S2,RHO) returns an array of random numbers chosen from
%   the Beckmann distribution
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
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BECKRND(U1,U2,S1,S2,RHO,M,N,...)   or
%   R = BECKRND(U1,U2,S1,S2,RHO,[M,N,...])
%   returns an M-by-N-by-... array.

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 5
    error('beckrnd:TooFewInputs','Requires at least five input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

mu = [u1 u2];
Sigma = [s1.^2 rho.*s1.*s2; rho.*s1.*s2 s2.^2];

sz=cell2mat(varargin);

num=prod(sz(:));


if isscalar(sz)
    num=num.^2;
    sz=[sz sz];
end

bvn = mvnrnd(mu, Sigma, num);
r=sqrt(bvn(:,1).^2+bvn(:,2).^2);
r=reshape(r,sz);

end