function [m,v]=gum1stat(a,b)
%GUM1STAT Mean and variance of the Type-1 Gumbel Distribution
%   [M,V]=GUM1STAT(A,B) returns the mean and variance of the Type-1 Gumbel
%   Distribution with parameters A and shape parameter B.
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('gum1stat:TooFewInputs',...
          'Requires two input arguments.'); 
end


%Initialize M to NaN.
if isa(a,'single') || isa(b,'single') 
    m=zeros(size(a),'single');
else
    m=zeros(size(a));
end
v=m;


try
    eg=0.57721566490153286060651209008240243104215933593992;
    m=(eg+log(b))./a;
    m2=((6*eg^2)+(pi^2)+((6*log(b)).*(2*eg+log(b))))./(6*(a.^2));
    v=m2-(m.^2);
catch
    error('gum1stat:InputSizeMismatch');
end



end