function y = daviscdf(x,b,n,u)
%NOTE: NEEDS ZETA FUNCTION, SYMBOLIC TOOLBOX
%DAVISCDF Davis probability density function
%   Y = DAVISCDF(X,B,N,U) returns the Davis probability density function
%   with with scale parameter B, shape parameter N, and location 
%   parameter U, at the values in X.
%
%   Restrictions:
%   B: any positive real number
%   N: great than 1
%   U: any non-negative real number
%   Valid for X>U

%   The size of Y is the common size of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%

%   Mike Sheppard
%   Last Modified 9-May-2011

if nargin < 4
    error('daviscdf:TooFewInputs',...
          'Requires at least four input argument.'); 
end

[errorcode x b n u] = distchck(4,x,b,n,u);

if errorcode > 0
    error('daviscdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(b,'single') || isa(n,'single') || isa(u,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


% Return NaN if not valid argument
y( b<=0 | n<=1 | u<0 | x<u ) = NaN;

k=find(b>0 & n>1 & u>=0 & x>u);
if any(k),
    xk = x(k);
    bk = b(k);
    nk = n(k);
    uk = u(k);
    sum2=zeros(size(xk));
    %No closed form, CDF expressed as Taylor series of integral with
    %substitution of variables
    %Let z=b/(t-u)   (t dummy variable for x)
    %Then CDF = 1 - (INTEGRAL / (Gamma[n]*Zeta[n]))
    %Where INTEGRAL = Integral([z^(n-1) / (e^z-1)] , dz, 0, b/(x-u)]
    %And can be expressed as sum of taylor series as
    %Sum( (B_k / k!*(n+k-1)) * z^(n+k-1) , k, 0 , Inf]
    %where B_k is the k'th Bernoulli number  
    zk=bk./(xk-uk);
    for j=0:50
        term=mfun('bernoulli',j).*(zk.^(j+nk-1))./(gamma(j+1).*(j+nk-1));
        %term=bern(j).*(zk.^(j+nk-1))./(gamma(j+1).*(j+nk-1));
        sum2=sum2+term;
    end
    y(k)=1-((sum2)./(gamma(nk).*zeta(nk)));

end


end