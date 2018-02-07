function x = triinv(p,xmin,xmax,c)
%TRIINV Inverse Cumulative Triangle Distribution probability density
%   X = TRICDF(P,XMIN,XMAX,C) returns the probability density of the
%   Triangle Distribution with parameters XMIN, XMAX and with a maximum
%   at C. 
%
%   If C is not given it is assumed to be symmetric with C at the midpoint.
%   That is, C=(XMIN+XMAX)/2

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 3
    error('triinv:TooFewInputs',...
          'Requires at least three input arguments.');
end

if nargin==3
    c=(xmin+xmax)./2;
end

[errorcode, p,xmin,xmax,c] = distchck(4,p,xmin,xmax,c);

if errorcode > 0
    error('triinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


if any((c<xmin)|(c>xmax))
    error('tricdf:ValidInputs',...
          'C must be betwen Xmin and Xmax.');
end

% Initialize P to NaN.
if isa(p,'single') || isa(xmin,'single') || isa(xmax,'single') || isa(c,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

cutoff=(c-xmin)./(xmax-xmin);
k1=(p>=0)&(p<=cutoff);
k2=(p>cutoff)&(p<=1);

if any(k1)
    temp=p(k1).*(c(k1)-xmin(k1)).*(xmax(k1)-xmin(k1));
    x(k1)=xmin(k1)+sqrt(temp);
end

if any(k2)
    temp=(1-p(k2)).*(-c(k2)+xmax(k2)).*(xmax(k2)-xmin(k2));
    x(k2)=xmax(k2)-sqrt(temp);   
end


end