function y = irwinhallpdf(x,n,a,b)
%IRWINHALLPDF Irwin-Hall probability density function
%   Y = IRWINHALLPDF(X,N,A,B) returns the probability density function of
%   the Irwin-Hall Distribution of the sum of N independent random
%   variables uniformly distributed from A to B.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Distribution: Continuous, bounded, (A*N,B*N)
%   Restrictions:
%      A <= (X\N) <= B
%      N >= 1        (integer)
%
%   NOTE: The Irwin-Hall Distribution is also known as the
%   Uniform Sum Distribution. For the mean of N uniformly distributed
%   random variables, use BATESPDF. Loss of accuracy occurs when N is
%   greater than 25.
%
%   See also IRWINHALLCDF, IRWINHALLINV, IRWINHALLSTAT, IRWINHALLFIT,
%            IRWINHALLLIKE, IRWINHALLRND, IRWINHALLSF, IRWINHALLHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin < 2
    error('irwinhallpdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin == 2
    a=0; b=1;
elseif nargin==3
    b=1;
end


[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('irwinhallpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>=1) & (n==round(n));
okvar = (x>=a.*n) & (x<=b.*n);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    a=a(ok); b=b(ok); n=n(ok); x=x(ok);
    
    
    sz=size(x);                          %size to reshape afterwards
    x=x(:); n=n(:); aa=a(:); bb=b(:);    %vectorize all
    aa(bb<aa)=NaN; bb(bb<aa)=NaN;        %Restrictions
    x=(x-n.*aa)./(bb-aa);                %Transform to standard sums of U[0,1]
    
    
    %Make everything matrix
    km=repmat(0:max(n),length(x),1);   %Matrix k
    xm=repmat(x,1,size(km,2));         %Matrix X
    nm=repmat(n,1,size(km,2));         %Matrix n
    
    %Pre-allocate
    M=zeros(size(km));  %pre-allocate
    nck=exp(gammaln(nm+1)-gammaln(km+1)-gammaln(nm-km+1)); %matrix nchoosek
    M=((-1).^km).*nck.*((xm-km).^(nm-1)).*sign(xm-km);
    %For loop, different number of sums per row
    sM=zeros(size(M,1),1);
    %Hopefully can be vectorized later
    if length(unique(n))==1
        sM=sum(M,2);
    else
        for i=1:length(x)
            sM(i)=sum(M(i,1:1+n(i)));
        end
    end
    sM=sM(:);
    yt=sM./(2*gamma(n));
    %Replace any n=1 with unifpdf
    k1=(n==1);
    if any(k1)
        yt(k1)=1; %Going to divide by (b-a) later on
    end
    
    
    yt=reshape(yt(:),sz);
    y(ok)=yt./(b-a);
    
end

%Round-off
y(y<0)=0;


end