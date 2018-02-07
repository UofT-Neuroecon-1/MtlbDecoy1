function [m,v] = poissbinostat(p)
%POISSBINOSTAT Mean and variance for Poisson binomial distribution
%   Y = POISSBINOSTAT(P) returns the mean and variance of the Poisson
%   Binomial distribution of the sum of independent Bernoulli trials given
%   by the probability vector P
%
%   P is a vector of length N, representing the N Bernoulli probabilities
%
%   If P(i)=p for all i=1,...,N then simplifies to Binomial Distribution
%   POISSBINOSTAT(P)=BINOSTAT(N,p)
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      x>=0 (integer)
%      0<=P(i)<=1, for all i
%
%   The size of Y is the size of the input variable X.
%

% CHECK DISTRIBUTION  (sum=1)

%   Mike Sheppard
%   Last Modified 14-Jun-2011


if nargin < 1
    error('poissbinostat:TooFewInputs',...
          'Requires one input argument.');
end


%P to row vector, for calculation purposes
p=(p(:))';


%Return error if any P is out of bounds
if any(p<0) | any(p>1)
    error('poissbinostat:Input',...
          'Input vector P must be between 0 and 1.');
end


% Initialize m to zero.
if isa(p,'single')
    m=zeros(size(p),'single');
else
    m=zeros(size(p));
end
v=m;

m=sum(p);
v=sum(p.*(1-p));

end