function x = radinv(p)
%RADINV Inverse of the Rademacher cumulative distribution function
%   x = RADINV(p) returns the inverse of the Rademacher cumulative
%   distribution function
%
%   x(p< .5) = -1
%   x(p>=.5) = +1
%
%   Type: Discrete, Bounded
%   Restrictions:
%   0<=p<=1
%
%   The size of x is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Jun-2011


if nargin < 1
    error('radinv:TooFewInputs',...
          'Requires one input argument.');
end

% Initialize X to zero.
if isa(p,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

x(p<0 | p>1)=NaN;
x(p>0 & p<=.5)=-1;
x(p>.5 & p<=1)=+1;

end
