function p = tricdf(x,xmin,xmax,c)
%TRICDF Cumulative Triangle Distribution probability density
%   P = TRICDF(X,XMIN,XMAX,C) returns the probability density of the
%   Triangle Distribution with parameters XMIN, XMAX and with a maximum
%   at C. 
%
%   If C is not given it is assumed to be symmetric with C at the midpoint.
%   That is, C=(XMIN+XMAX)/2

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 3
    error('tricdf:TooFewInputs',...
          'Requires at least three input arguments.');
end

if nargin==3
    c=(xmin+xmax)./2;
end

[errorcode, x,xmin,xmax,c] = distchck(4,x,xmin,xmax,c);

if errorcode > 0
    error('tricdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


if any((c<xmin)|(c>xmax))
    error('tricdf:ValidInputs',...
          'C must be betwen Xmin and Xmax.');
end

% Initialize P to NaN.
if isa(x,'single') || isa(xmin,'single') || isa(xmax,'single') || isa(c,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end


k1=(x>=xmin)&(x<=c);
k2=(x>c)&(x<=xmax);

if any(k1)
    num=(-xmin(k1)+x(k1)).^2;
    den=(c(k1)-xmin(k1)).*(xmax(k1)-xmin(k1));
    p(k1)=num./den;
end

if any(k2)
    num=(xmax(k2)-x(k2)).^2;
    den=(-c(k2)+xmax(k2)).*(xmax(k2)-xmin(k2));
    p(k2)=1-(num./den);    
end

%Out of bounds
p(x<xmin)=0; p(x>xmax)=1; 

end
