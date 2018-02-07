function x = gum1inv(p,a,b)
%GUM1INV Inverse of Type-1 Gumbel distribution
%   X = GUM1INV(P,A,B) returns the inverse of Type-1 Gumbel distribution
%   with parameters A and shape parameter B
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 3
    error('gum1inv:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('gum1inv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

k1=(p>0)&(p<1);
if any(k1)
    temp=-log(p(k1))./b(k1);
    x(k1)=-log(temp)./a(k1);
end

%End cases
x(p==0)=0;
x(p==1)=Inf;


end