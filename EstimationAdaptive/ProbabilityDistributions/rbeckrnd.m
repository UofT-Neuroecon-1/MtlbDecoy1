function r = rbeckrnd(s1,s2,varargin)
%RBECKRND Random arrays from the Reduced Beckmann distribution
%   R = RBECKRND(S1,S2) returns an array of random numbers chosen from
%   the Reduced Beckmann distribution
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
%   R = RBECKRND(S1,S2,M,N,...) or R = RBECKRND(S1,S2,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 3-Dec-2011


if nargin < 2
    error('rbeckrnd:TooFewInputs','Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Use definition of Beckmann distribution
sigma = [s1.^2 s2.^2];
sz=cell2mat(varargin);

num=prod(sz(:)); %number of random numbers

if isscalar(sz) %If scalar then represents a square array
    num=num.^2;
    sz=[sz sz];
end

bvn = mvnrnd([0 0], sigma, num);
r=sqrt(bvn(:,1).^2+bvn(:,2).^2); %Compute distance
r=reshape(r,sz); %Reshape to proper dimensions


end