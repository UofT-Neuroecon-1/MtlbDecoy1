function [m,v] = sgomstat(b,n)
%SGOSTAT Mean and Variance
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


if nargin < 2
    error('sgomstat:TooFewInputs',...
          'Requires two input arguments.');
end


[errorcode,b,n] = distchck(2,b,n);

if errorcode > 0
    error('sgomstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(b,'single') || isa(n,'single')
    m=zeros(size(p),'single');
else
    m=zeros(size(p));
end
v=m;

%Out of bounds
b(b<0)=NaN; n(n<0)=NaN;
%Check to see if 'lambertw' exists; if not do Newton's Method

chi=mfun('Chi', n);
shi=mfun('Shi', n);
eulergamma=0.57721566490153286060651209008240243104215933593992;
m=(1+eulergamma.*n-cosh(n)+sinh(n)+n.*(-chi+log(n)+shi(n)))./(b.*n);


end