function p = gum1cdf(x,a,b)
%GUM1CDF Cumulative Type-1 Gumbel distribution
%   P = GUM1CDF(X,A,B) returns the cumulative Type-2 Gumbel density
%   function with parameters A and Shape Parameter B
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 3
    error('gum1cdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('gum1cdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end


term1=exp(-a.*x);
term2=-b.*term1;
p=exp(term2);

end