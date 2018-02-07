function p = pertcdf(x,xmin,xmax,c,lambda)
%PERTCDF PERT cumulative distribution function
%   P = PERTCDF(X,XMIN,XMAX,C,LAMBDA) returns the PERT cumulative 
%   distribution fun with range XMIN to XMAX and maximum at C with
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
%      XMIN<=X<=XMAX  Bounded
%      XMIN<C<XMAX    With maximum within the bounds
%      LMABDA>0       non-negative scale


%   Mike Sheppard
%   Last Modified 1-Jul-2011

if nargin < 1
    error('pertcdf:TooFewInputs',...
          'Requires at least one input.');
elseif nargin==1
    xmin=0; xmax=1; c=.5; lambda=1;
elseif nargin==2
    error('pertcdf:TooFewInputs',...
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
    error('pertcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Transformation of the Beta Distribution
xdiff=xmax-xmin;
trans_x=(x-xmin)./xdiff;
a=1+lambda.*((c-xmin)./xdiff);
b=1+lambda.*((xmax-c)./xdiff);
p=betacdf(trans_x,a,b);

end