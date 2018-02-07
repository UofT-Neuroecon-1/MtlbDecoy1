function [m,v]=vonmstat(u,k,ind)
%VONMSTAT Mean and variance of the Von-Mises distribution
%   [M,V]=VONMSTAT(U,K) returns the mean of the Von-Mises distribution with
%   mean U and concentration K.
%
%   If IND=0 the expected value is M=U, by definition. Variance is not
%   computed
%   If IND=1 the circular mean and variance are computed
%

%   The size of the output is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 2
    error('vonmstat:TooFewInputs',...
          'Requires at least two input arguments.');
end

if nargin==2
    ind=0; %Regular expected value
end

[errorcode, u,k,ind] = distchck(3,u,k,ind);

if errorcode > 0
    error('vonmstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Initialize x to NaN
if isa(u,'single') || isa(k,'single')
   m = zeros(size(u),'single');
else
   m = zeros(size(u));
end
v=m;

rm=(ind==0);
cm=(ind==1);
if any(rm)
    m(rm)=u(rm); %by definition
end
if any(cm)
    m(cm)=besseli(1,k(cm))./besseli(0,k(cm));
    v(cm)=1-m;
end


end