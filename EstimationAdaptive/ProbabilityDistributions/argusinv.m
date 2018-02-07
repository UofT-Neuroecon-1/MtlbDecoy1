function x = argusinv(p,chi,c)
%ARGUSINV Inverse of the ARGUS cumulative distribution function
%   X = ARGUSINV(P,CHI,C) returns the inverse cumulative distribution
%   function of the ARGUS Distribution with curvature parameter CHI,
%   and cut-off parameter C, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1.
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%      CHI, C > 0
%
%   NOTE: ARGUSINV uses Newton's method to converge to the solution.
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSSTAT, ARGUSFIT, 
%            ARGUSLIKE, ARGUSRND, ARGUSSF, ARGUSHAZ
%

%   Mike Sheppard
%   Last Modified: 14-Dec-2011


if nargin < 2
    error('argusinv:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin == 2
    c=1;
end

[errorcode p chi c] = distchck(3,p,chi,c);

if errorcode > 0
    error('argusinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize x to zero
if isa(p,'single') || isa(chi,'single') || isa(c,'single')
    x=zeros(size(p),'single');
    seps=sqrt(eps('single'));
else
    x=zeros(size(p));
    seps=sqrt(eps);
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<chi & chi<Inf);
okvar = (0<p & p<1);
ok = (okparam & okvar);
x(~ok)=NaN;
%The inverse of cdf of 0 is 0, and the inverse cdf of 1 is c
x(okparam & p==0)=0;
x(okparam & p==1)=c(okparam & p==1);


if all(~ok), return; end

p=p(ok); chi=chi(ok); c=c(ok);


%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;


% Use the mean as a starting guess
[xok,ignore]=argusstat(chi,c);
if isa(p,'single')
    xok=single(xok);
end

%Move starting values away from the boundaries
xok(xok==0)=seps;
xok(xok==1)=1-seps;

h=ones(size(p));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(xok)) && count<count_limit),
    count=count+1;
    h=(arguscdf(xok,chi,c)-p)./arguspdf(xok,chi,c);
    xnew=xok-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = (xnew<=0);
    klarge = (xnew>=c);
    if any(ksmall) || any(klarge)
        xnew(ksmall)=xok(ksmall)/10;
        xnew(klarge)=c(klarge)-(c(klarge)-(xok(klarge)))/10;
    end
    xok=xnew;
end

%Return the converged value(s).
x(ok)=xok;

if count==count_limit
    fprintf('\nWarning: ARGUSINV did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end



end