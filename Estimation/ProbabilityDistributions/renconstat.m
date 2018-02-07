function [m,v] = renconstat(n)
%RENCONSTAT Mean and Variance for the Rencontres Distribution
%   [M,V] = RENCONSTAT(n) returns the mean and variance of the Recontres
%   Distribution of having X fixed points in a uniformly distributed random
%   permutation of { 1, ..., N }.
%
%   Type: Discrete, Bounded
%   Restrictions:
%     n>=1    (integer)
%
%   The size of the outputs is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011


if nargin < 1
    error('renconstat:TooFewInputs',...
          'Requires one input argument.');
end


% Initialize Y to zero.
if isa(n,'single')
    m=zeros(size(n),'single');
else
    m=zeros(size(n));
end
v=m;

m(n<1 | n~=round(n))=NaN;
v(n<1 | n~=round(n))=NaN;


k=(n>=1 & n==round(n));
if any(k)
   m(k)=1;
   v(k)=1; %Second moment is 2, v=2-1^2
end


end