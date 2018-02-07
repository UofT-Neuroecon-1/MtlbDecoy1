function p = slashcdf(x)
%SLASHCDF Slash Distribution
%   P = SLASHCDF(X) returns the probability density of

%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.  
%

%   Mike Sheppard
%   Last Modified 6-Jun-2011


if nargin < 1
    error('slashcdf:TooFewInputs',...
          'Requires one input arguments.'); 
end


% Initialize Y to zero.
if isa(x,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

term1=normcdf(x);
term2=(normpdf(0)-normpdf(x))./(x);
p=term1-term2;
%Removable discontinuity
p(x==0)=1/2;


end
