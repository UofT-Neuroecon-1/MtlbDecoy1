function x = davisinv(p,b,n,u)
%DAVISINV Inverse of the Davis cumulative distribution function (CDF)
%   X = DAVISINV(P,B,N,U) returns the inverse of the Davis Distribution
%   with scale parameter B, shape parameter N, and locaiton parameter U, at
%   the values in P.
%
%   Restrictions:
%   B: any positive real number
%   N: great than 1
%   U: any non-negative real number
%   Valid for X>U
%
%   The size of Y is the common size of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   NOTE: DAVISINV uses Newton's method to converge to the solution.
%   IT IS SLOW TO CONVERGE

%   Mike Sheppard
%   Last Modified 9-May-2011


if nargin < 4
    error('davisinv:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode p b n u] = distchck(4,p,b,n,u);

if errorcode > 0
    error('davisinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Initialize x to zero
if isa(p,'single') || isa(b,'single') || isa(n,'single') || isa(u,'single')
    x=zeros(size(p),'single');
    seps=sqrt(eps('single'));
else
    x=zeros(size(p));
    seps=sqrt(eps);
end

% Return NaN for out of range parameters.
x( p<0 | p>1 | b<=0 | n<=1 | u<0)=NaN;

%The inverse of cdf of 0 is u, and the inverse cdf of 1 is Inf
k0=find(p==0 & b>0 & n>1 & u>=0);
if any(k0)
    x(k0)=u(k0);
end
k1=find(p==1 & b>0 & n>1 & u>=0);
if any(k1);
    x(k1)=Inf;
end


%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;

k=find( p>0 & p<1 & b>0 & n>1 & u>=0);
if isempty(k)
    return;
end
pk=p(k);

% Use the mean as a starting guess
xk = davisstat(b(k),n(k),u(k));
if isa(p,'single')
    xk=single(xk);
end

%Move starting values away from the boundaries
xk(xk==u(k))=u(xk==u(k))+seps;
h=ones(size(pk));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(xk)) && max(abs(j))>crit && count<count_limit),
    count=count+1;
    h=(daviscdf(xk,b(k),n(k),u(k))-pk)./davispdf(xk,b(k),n(k),u(k));
    xnew=xk-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = find(xnew<=u(k));
    if any(ksmall)
        xnew(ksmall)=u(ksmall)+((u(ksmall)-xk(ksmall))/10);
    end
    xk=xnew;
end

%Return the converged value(s).
x(k)=xk;

if count==count_limit
    fprintf('\nWarning: DAVISINV did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end


%Reshape to original input
x=reshape(x,size(p));

end