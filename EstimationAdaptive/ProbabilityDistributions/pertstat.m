function [m,v]=pertstat(xmin,xmax,c,lambda)
%PERTSTAT Mean and variance for the PERT Distribution
%   [M,V] = PERTSTAT(XMIN,XMAX,C,LAMBDA) returns the mean and variance
%   of the PERT Distribution with range XMIN to XMAX and maximum at C with
%   shape parameter LAMBDA.
%
%   If LAMBDA is not provided it will default to 1.
%   If C and LAMBDA are not provided a symmetric distribution with lambda=1
%   will be computed
%   If XMIN, XMAX, C, and LAMBDA are not provided it will default to a
%   symmetric distribution from [0,1] with lambda=1
%
%   The PERT distribution is a reparametrized Beta Distribution.
%
%   Type: Continuous, bounded
%   Restrictions:
%      XMIN<C<XMAX    With maximum within the bounds
%      LMABDA>0       non-negative scale


%   Mike Sheppard
%   Last Modified 1-Jul-2011




if nargin==0
    xmin=0; xmax=1; c=.5; lambda=1;
elseif nargin==1
    error('pertstat:TooFewInputs',...
          'Requires at least two inputs.');
elseif nargin==2
    c=(xmax+xmin)./2; lambda=4;
elseif nargin==3
    lambda=4;
else
    %good
end


[errorcode xmin xmax c lambda] = distchck(4,xmin,xmax,c,lambda);

if errorcode > 0
    error('pertstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(xmin,'single') || isa(xmax,'single') || isa(c,'single') || isa(lambda,'single')
    m=zeros(size(xmin),'single');
else
    m=zeros(size(xmin));
end
v=m;

%Mean
m=(xmax+xmin+c.*lambda)./(2+lambda);

%Variance
num=(xmax-xmin-c.*lambda+lambda.*xmax).*(xmax+c.*lambda-xmin.*(1+lambda));
den=(2+lambda).^2.*(3+lambda);
v=num./den;

end
