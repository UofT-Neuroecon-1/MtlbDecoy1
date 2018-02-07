function [m,v]=gum2stat(a,b)
%GUM2STAT Mean and variance of the Type-2 Gumbel Distribution
%   [M,V]=GUM2STAT(A,B) returns the mean and variance of the Type-2 Gumbel
%   Distribution with parameters A and shape parameter B.
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('gum2stat:TooFewInputs',...
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
    m=((b.^(1./a)).*gamma(-1./a))./a;
    m2=-(b.^(2./a)).*gamma((-2+a)./a);
    v=m2-(m.^2);
catch
    error('gum2stat:InputSizeMismatch');
end



end