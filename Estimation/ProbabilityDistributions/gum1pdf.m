function y = gum1pdf(x,a,b)
%GUM1PDF Probability density function for the Type-1 Gumbel Distribution
%   Y = GUM1PDF(X,A,B) returns the probability density function of the
%   Type-1 Gumbel density function with parameters A and Shape Parameter B,
%   evaluated at the values in X.
%
%   Type: Continuous, unbounded
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 3
    error('gum1pdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('gum1pdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


term1=-(b.*exp(-a.*x)+(a.*x));
y=a.*b.*exp(term1);


end