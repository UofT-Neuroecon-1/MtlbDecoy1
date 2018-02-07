function y = tripdf(x,xmin,xmax,c)
%TRIPDF Triangle probability density function
%   Y = TRIPDF(X,XMIN,XMAX,C) returns the probability density of the
%   Triangle Distribution with parameters XMIN, XMAX and with a maximum
%   at C, at the values in X.
%
%   TRIPDF(X,XMIN,XMAX) assumes symmetric distribution at the midpoint.
%
%   Type: Continuous, bounded
%   Restrictions:
%        XMIN <= (X, C) <= XMAX
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%


%   Mike Sheppard
%   Last Modified 31-Jul-2011


if nargin < 3
    error('tripdf:TooFewInputs',...
          'Requires at least three input arguments.');
end

if nargin==3
    c=(xmin+xmax)./2;
end

[errorcode, x,xmin,xmax,c] = distchck(4,x,xmin,xmax,c);

if errorcode > 0
    error('tripdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


if any((c<xmin)|(c>xmax))
    error('tripdf:ValidInputs',...
          'C must be betwen Xmin and Xmax.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(xmin,'single') || isa(xmax,'single') || isa(c,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


k1=(x>=xmin)&(x<=c);
k2=(x>c)&(x<=xmax);

if any(k1)
    num=2.*(-xmin(k1)+x(k1));
    den=(c(k1)-xmin(k1)).*(xmax(k1)-xmin(k1));
    y(k1)=num./den;
end

if any(k2)
    num=2.*(xmax(k2)-x(k2));
    den=(-c(k2)+xmax(k2)).*(xmax(k2)-xmin(k2));
    y(k2)=num./den;    
end

end

