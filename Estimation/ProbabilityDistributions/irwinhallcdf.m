function p = irwinhallcdf(x,n,a,b)
%IRWINHALLCDF Irwin-Hall cumulative distribution function
%   P = IRWINHALLCDF(X,N,A,B) returns the cumulative distribution function
%   of the Irwin-Hall Distribution of the sum of N independent random
%   variables uniformly distributed from A to B.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Type: Continuous, bounded, (A*N,B*N)
%   Restrictions:
%        A <= (X\N) <= B
%        N >= 1        (integer)
%
%   Note: The Irwin-Hall Distribution is also known as the
%   Uniform Sum Distribution. For the mean of N uniformly distributed
%   random variables, use BATESPDF. Loss of accuracy occurs when N is
%   greater than 25.
%
%   See also IRWINHALLPDF, IRWINHALLINV, IRWINHALLSTAT, IRWINHALLFIT,
%            IRWINHALLLIKE, IRWINHALLRND, IRWINHALLSF, IRWINHALLHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012

%ERROR:
% p = irwinhallcdf(linspace(-3,6),3,-1,1); plot(p)
 
if nargin < 2
    error('irwinhallcdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin == 2
    a=0; b=1;
elseif nargin==3
    b=1;
end


[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('irwinhallcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>=1) & (n==round(n));
okvar = (x>=a.*n) & (x<=b.*n);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<a.*n)=0;
p(okparam & x>b.*n)=1;


if any(ok)
    a=a(ok); b=b(ok); n=n(ok); x=x(ok);
    
    
    sz=size(x);                          %size to reshape afterwards
    x=x(:); n=n(:); aa=a(:); bb=b(:);    %vectorize all
    x=(x-n.*aa)./(bb-aa);                %Transform to standard sums of U[0,1]
    
    
    %Make everything matrix
    km=repmat(0:max(n),length(x),1);   %Matrix k
    xm=repmat(x,1,size(km,2));         %Matrix X
    nm=repmat(n,1,size(km,2));         %Matrix n
    
    %Pre-allocate
    M=zeros(size(km));  %pre-allocate
    nck=exp(gammaln(nm+1)-gammaln(km+1)-gammaln(nm-km+1)); %matrix nchoosek
    %M=((-1).^km).*nck.*((xm-km).^(nm-1)).*sign(xm-km);  %pdf
    %Closed form expression for the CDF involves Heaviside functions
    %Which are the integral of the last two terms of the pdf as a function of x;
    %after computing Heaviside terms, multiple by first two terms and add.
    temp_hsf=(2.*((-km+xm).^nm).*hsf(-km+xm).*(hsf(xm)+hsf(-km,-xm)))-...
        (2.*((-km).^nm).*hsf(-km).*(hsf(-xm)+hsf(-xm,-km+xm)));
    temp_inside=(1./nm).*(((-1).^nm.*((km.^nm)-(km-xm).^nm))+temp_hsf);
    M=((-1).^km).*nck.*temp_inside;  %cdf
    
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
    pt=sM./(2*gamma(n));
    %Replace any n=1 with unifpdf
    k1=(n==1);
    if any(k1)
        pt(k1)=x(k1); %Going to divide by (b-a) later on
    end
    
    
    pt=reshape(pt(:),sz);
    p(ok)=pt./(b-a);
    
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end



function c=hsf(varargin)
%Multidimensional Heaviside function, which is 1 only if none of
%the inputs are not positive
for i1=1:length(varargin)
    M1(:,:,i1)=varargin{i1};
end
c=all(M1>0,3);
end