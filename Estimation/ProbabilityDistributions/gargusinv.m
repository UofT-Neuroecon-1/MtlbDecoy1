function x = gargusinv(y,p,chi,c)
%GARGUSINV Inverse of the Generalized ARGUS cdf
%   X = GARGUSINV(Y,P,CHI,C) returns the inverse cdf of the Generalized 
%   ARGUS Distribution with power P, curvature CHI, and cut-off C, 
%   at the values in Y.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1.
%
%   Type: Continuous, bounded, (0,C)
%   Restrictions:
%        P > -1
%        CHI, C > 0
%
%   Note: GARGUSINV uses Newton's method to converge to the solution.
%
%   See also GARGUSPDF, GARGUSCDF, GARGUSSTAT, GARGUSFIT, 
%            GARGUSLIKE, GARGUSRND, GARGUSSF, GARGUSHAZ
%

%   Mike Sheppard
%   Last Modified: 20-Dec-2011


if nargin < 3
    error('gargusinv:TooFewInputs',...
        'Requires at least three input arguments.');
end

if nargin == 3
    c=1;
end

[errorcode y p chi c] = distchck(4,y,p,chi,c);

if errorcode > 0
    error('argusinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize x to zero
if isa(y,'single') || isa(p,'single') || isa(chi,'single') || isa(c,'single')
    x=zeros(size(y),'single');
    seps=sqrt(eps('single'));
else
    x=zeros(size(y));
    seps=sqrt(eps);
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-1<p & p<Inf) & (0<c & c<Inf) & (0<chi & chi<Inf);
k = (okparam & (0 < y & y < 1));
x(~k)=NaN;
%The inverse of cdf of 0 is 0, and the inverse cdf of 1 is c
x(okparam & y==0)=0;
x(okparam & y==1)=c(okparam & y==1);


if isempty(k)
    return;
end
yk=y(k);

%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;


% Use the mean as a starting guess
[xk,ignore]=gargusstat(p(k),chi(k),c(k));
if isa(p,'single')
    xk=single(xk);
end

%Move starting values away from the boundaries
xk(xk==0)=seps;
xk(xk==1)=1-seps;

h=ones(size(yk));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(xk)) && count<count_limit),
    count=count+1;
    h=(garguscdf(xk,p(k),chi(k),c(k))-yk)./garguspdf(xk,p(k),chi(k),c(k));
    xnew=xk-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = (xnew<=0);
    ck=c(k);
    klarge = (xnew>=ck);
    if any(ksmall) || any(klarge)
        xnew(ksmall)=xk(ksmall)/10;
        xnew(klarge)=ck(klarge)-(ck(klarge)-(xk(klarge)))/10;
    end
    xk=xnew;
end

%Return the converged value(s).
x(k)=xk;

if count==count_limit
    fprintf('\nWarning: GARGUSINV did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end


%Reshape to original input
x=reshape(x,size(y));



end