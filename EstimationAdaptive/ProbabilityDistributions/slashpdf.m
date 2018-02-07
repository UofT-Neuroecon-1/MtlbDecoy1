function y = slashpdf(x)
%SLASHPDF Slash Distribution probability density function
%   Y = SLASHPDF(X) returns the probability density function of the Slash
%   Distribution at the values in X
%
%   Type: Continuous, unbounded
%   Restrictions:
%        X real
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs. 
%

%   Mike Sheppard
%   Last Modified 5-Jul-2011


if nargin < 1
    error('slashpdf:TooFewInputs',...
        'Requires one input argument.');
end


% Initialize Y to zero.
if isa(x,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

k=(x~=0);
if any(k)
    y(k)=(normpdf(0)-normpdf(x(k)))./(x(k).^2);
end

%Removable discontinuity
y(x==0)=normpdf(0)/2;


end
