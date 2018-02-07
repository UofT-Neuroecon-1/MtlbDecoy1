function y = gum2pdf(x,a,b)
%GUM2PDF Probability density function for the Type-2 Gumbel Distribution
%   Y = GUM2PDF(X,A,B) returns the probability density function of the
%   Type-2 Gumbel density function with parameters A and Shape Parameter B,
%   evaluated at the values in X.
%
%   Type: Continuous, unbounded
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 3
    error('gum2pdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('gum2pdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


term1=a.*b.*(x.^(-a-1));
term2=-b.*(x.^(-a));
y=term1.*exp(term2);


end