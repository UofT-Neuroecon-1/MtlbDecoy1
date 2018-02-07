function x = irwinhallinv(p,n,a,b)
%IRWINHALLINV Inverse of the Irwin-Hall cumulative distribution function (cdf)
%   X = IRWINHALLINV(P,N,A,B) returns the inverse of the Irwin-Hall
%   distribution of the sum of N independent and identically distributed
%   random variables uniformly distributed continuously from A to B
%
%   Default values for A=0 and B=1.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   IRWINHALLINV uses Newton's method to converge to the solution

%   Mike Sheppard
%   Last Modified: 21-Mar-2011

%   Algorithm based on 'betainv' revision 2.11.2.5 (2004/12/06)
%   Which uses the Newton's Method



if nargin < 2
    error('irwinhallinv:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin==2
    a=0;
    b=1;
end


[errorcode p n a b] = distchck(4,p,n,a,b);

if errorcode > 0
    error('irwinhallinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize x to zero
if isa(p,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
    seps=sqrt(eps('single'));
else
    x=zeros(size(p));
    seps=sqrt(eps);
end

% Return NaN for out of range parameters.
x( p<0 | p>1 | n<1 | n~=round(n))=NaN;

%The inverse of cdf of 0 is a*n, and the inverse cdf of 1 is b*n
k0=find(p==0 & n==round(n));
if any(k0)
    x(k0)=a(k0).*n(k0);
end
k1=find(p==1);
if any(k1);
    x(k1)=b(k1).*n(k1);
end


%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;

k=find( p>0 & p<1 & n==round(n));
if isempty(k)
    return;
end
pk=p(k);

% Use the mean as a starting guess
xk = n(k).*(a(k)+b(k))/2;
if isa(p,'single')
    xk=single(xk);
end

%Move starting values away from the boundaries
mink=n(k).*a(k);
maxk=n(k).*b(k);
xk(xk==mink)=seps;
xk(xk==maxk)=1-seps;

h=ones(size(pk));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(xk)) && max(abs(j))>crit && count<count_limit),
    count=count+1;
    h=(irwinhallcdf(xk,n(k),a(k),b(k))-pk)./irwinhallpdf(xk,n(k),a(k),b(k));
    xnew=xk-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = find(xnew<=mink);
    klarge = find(xnew>=maxk);
    if any(ksmall) || any(klarge)
        xnew(ksmall)=xk(ksmall)/10;
        xnew(klarge)=maxk-(maxk-(xk(klarge)))/10;
    end
    xk=xnew;
end

%Return the converged value(s).
x(k)=xk;

if count==count_limit
    fprintf('\nWarning: IRWINHALLINV did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end


%Reshape to original input
x=reshape(x,size(p));

end