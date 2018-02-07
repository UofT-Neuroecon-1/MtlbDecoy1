function y = pertpdf(x,xmin,xmax,c,lambda)
%PERTPDF PERT probability density function
%   Y = PERTPDF(X,XMIN,XMAX,C,LAMBDA) returns the probability density 
%   function of the PERT Distribution with range XMIN to XMAX, 
%   maximum at C, shape parameter LAMBDA, at the values in X.
%
%   PERT(X) returns the symmetric PERT distribution on (0,1)
%   PERT(X,XMIN,XMAX) returns the symmetric PERT distribution
%                     on (XMIN,XMAX)
%   PERT(X,XMIN,XMAX,C) is the same as PERT(X,XMIN,XMAX,C,1)
%
%   Type: Continuous, bounded
%   Restrictions:
%      XMIN <= X <= XMAX 
%      XMIN < C < XMAX
%      LAMBDA > 0
%
%   Note: The PERT distribution is a reparametrization
%   of the Beta Distribution.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011

if nargin < 1
    error('pertpdf:TooFewInputs',...
          'Requires at least one input.');
elseif nargin==1
    xmin=0; xmax=1; c=.5; lambda=1;
elseif nargin==2
    error('pertpdf:TooFewInputs',...
          'Requires at least three inputs.');
elseif nargin==3
    c=(xmax+xmin)./2; lambda=4;
elseif nargin==4
    lambda=4;
else
    %good
end


[errorcode x xmin xmax c lambda] = distchck(5,x,xmin,xmax,c,lambda);

if errorcode > 0
    error('pertpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN for out of range parameters.
x(x<xmin | x>xmax)=NaN; xmin(xmin>xmax)=NaN; xmax(xmax<xmin)=NaN;
c(c<xmin | c<xmax)=NaN; lambda(lambda<0)=NaN;


%Transformation of the Beta Distribution
xdiff=xmax-xmin;
trans_x=(x-xmin)./xdiff;
a=1+lambda.*((c-xmin)./xdiff);
b=1+lambda.*((xmax-c)./xdiff);
y=betapdf(trans_x,a,b)./xdiff;

%Round off
y(y<0)=0;

end