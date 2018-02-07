function p = bortancdf(x,alpha,n)
%BORTANCDF Borel-Tanner cumulative distribution function
%   P = BORTANCDF(X,A,N) returns the cumulative distribution function
%   of the Borel-Tanner Distribution with shape parameters ALPHA and N,
%   at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, semi-bounded, {N,...,Inf}
%   Restrictions:
%        N >= 1         (N integer)
%        0 < ALPHA < 1
%
%   See also BORTANPDF, BORTANINV, BORTANSTAT, BORTANFIT, 
%            BORTANLIKE, BORTANRND, BORTANSF, BORTANHAZ
%


%   Mike Sheppard
%   Last Modified 14-May-2012



if nargin~=3
    error('bortancdf:TooFewInputs','Requires three input arguments.');
end

[errorcode x alpha n] = distchck(3,x,alpha,n);

if errorcode > 0
    error('bortancdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(alpha,'single') || isa(n,'single') 
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<alpha & alpha<1) & (n>=1 & n==round(n));
x=floor(x); %Floor for CDF for discrete
okvar = (n <= x) & (x < Inf) & (x==round(x));
ok=(okparam & okvar);


scalarparam=isscalar(alpha);
%This method sums the PDF to the required amount.
if any(ok)
    x=x(ok); alpha=alpha(ok); n=n(ok);
    val=max(x(:));
    I=(0:val)';
    if scalarparam
        tmp=cumsum(bortanpdf(I,alpha(1),n(1)));
        p(ok)=tmp(floor(x+1));  %CDF is right continuous
    else
        n1=sum(ok); n2=numel(I);
        compare=repmat(I,1,n1);
        index=repmat(x,n2,1);
        alphabig=repmat(alpha,n2,1);
        nbig=repmat(n,n2,1);
        y0=bortanpdf(compare,alphabig,nbig);
        y0(compare>index)=0;
        p(ok)=sum(y0,1);
    end
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end