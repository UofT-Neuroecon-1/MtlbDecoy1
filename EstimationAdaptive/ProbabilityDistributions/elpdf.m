function y = elpdf(x,p,b)
%ELPDF Exponential-logarithmic density function
%   Y = ELPDF(X,P,B) returns the Exponential-logarithmic density function
%   with paramaters P (0<P<1) and B (B>0) at the values in X.
%
%   The size of Y is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%

%   Mike Sheppard
%   Last Modified 10-May-2011

if nargin < 3
    error('elpdf:TooFewInputs',...
          'Requires at least three input argument.'); 
end

[errorcode x p b] = distchck(3,x,p,b);

if errorcode > 0
    error('elpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to NaN.
if isa(x,'single') || isa(p,'single') || isa(b,'single')
    y=NaN(size(x),'single');
else
    y=NaN(size(x));
end


k1=find(x>0 & p>0 & p<1 & b>0);
if any(k1)
    xk = x(k1);
    pk = p(k1);
    bk = b(k1);
    num=bk.*(1-pk).*exp(-bk.*xk);
    den=1-(1-pk).*exp(-bk.*xk);
    y(k1)=(num./den).*(1./(-log(pk)));
end

end