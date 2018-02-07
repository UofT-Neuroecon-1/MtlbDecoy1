function x = sgominv(p,b,n)
%SGOMINV Inverse of the Shifted Gompertz cumulative distribution function
%   x = SGOMINV(p,b,n) returns the inverse of the shifted gompertz 
%   probability density
%   function with scale parameter B and shape parameter N at the values in
%   X.
%
%   P,B,N are all real numbers great than 0.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 10-Jun-2011


if nargin < 3
    error('sgominv:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,p,b,n] = distchck(3,p,b,n);

if errorcode > 0
    error('sgominv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(p,'single') || isa(b,'single') || isa(n,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

%Out of bounds
b(b<0)=NaN; n(n<0)=NaN;
%Check to see if 'lambertw' exists; if not do Newton's Method

term1=lambertw(n.*p.*exp(n));
term2=log(n./(n-term1));
x=term2./b;


end