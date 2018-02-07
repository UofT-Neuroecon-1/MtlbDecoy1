function y = ncxpdf(x,v,delta)
%NCXPDF Non-Central Chi probability density function
%   Y = NCXPDF(X,V,DELTA) returns the probability density function of the
%   Non-Central Chi Distribution with V degrees of freedom, non-centrality
%   parameter, DELTA, at the values in X.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     X, V, DELTA > 0
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 3
    error('ncxpdf:TooFewInputs','Requires three input arguments.'); 
end

[errorcode x v delta] = distchck(3,x,v,delta);

if errorcode > 0
    error('ncxpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(v,'single') || isa(delta,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

k = (x>=0 & v>0 & delta>0);
if any(k)
    xk=x(k);
    vk=v(k);
    dk=delta(k);
    bes = besseli((vk/2)-1,dk.*xk,0);
    y(k) = exp((-0.5*(xk.^2+dk.^2))+vk.*log(xk)+log(dk)-(vk/2).*log(dk.*xk)).*bes;
end

% Return NaN for out of range parameters.
y(x<0 | v<=0 | delta<=0)=NaN;

%Round off
y(y<0)=0;

end