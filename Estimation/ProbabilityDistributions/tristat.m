function [m,v]=tristat(xmin,xmax,c)
%TRISTAT Mean and variance of the Triangle Distribution
%   [M,V]=TRISTAT(XMIN,XMAX,C) returns the mean and variance of the
%   Triangle Distribution with parameters XMIN, XMAX and with a maximum
%   at C.
%
%   If C is not given it is assumed to be symmetric with C at the midpoint.
%   That is, C=(XMIN+XMAX)/2

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 2
    error('tristat:TooFewInputs',...
          'Requires at least two input arguments.');
end

if nargin==2
    c=(xmin+xmax)./2;
end

[errorcode, xmin,xmax,c] = distchck(3,xmin,xmax,c);

if errorcode > 0
    error('tristat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


if any((c<xmin)|(c>xmax))
    error('tricdf:ValidInputs',...
          'C must be betwen Xmin and Xmax.');
end

% Initialize M and V to NaN.
if isa(xmin,'single') || isa(xmax,'single') || isa(c,'single')
    m=NaN(size(c),'single');
else
    m=NaN(size(c));
end
v=m;

m=(c+xmax+xmin)/3; %If symmetric reduces to (xmin+xmax)/2
v=(c.^2-(c.*xmax)+(xmax.^2)-(c.*xmin)-(xmax.*xmin)+(xmin.^2))/18;


end