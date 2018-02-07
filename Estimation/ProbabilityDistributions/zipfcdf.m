function y = zipfcdf(x,p,N)
%ZIPFCDF Cumulative Zipf probability density function
%   Y = ZIPFCDF returns the cumulative Zipf density function with
%   with parameter P and range N.
%
%   If N is not given default is Infinity
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 24-May-2011


if nargin < 2
    error('zipfcdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end

if nargin==2
    N=Inf;
end

[errorcode x p N] = distchck(3,x,p,N);

if errorcode > 0
    error('zipfcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(p,'single') || isa(N,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end



%Valid cases for finite case
k1=isfinite(N)&(p>0)&(rem(x,1)==0)&(x<=N)&(x>=1);
if any(k1)
    %For each k1, create harmonic number vector
    %(Should be able to do it without a for-loop)
    kin=find(k1);
    for i=1:length(kin)
    hmvn(i)=sum(1./((1:floor(x(kin(i)))).^(1+p(kin(i)))));   
    hmvd(i)=sum(1./((1:N(kin(i))).^(1+p(kin(i)))));
    end
    y(k1)=hmvn./hmvd;
end

%Valid cases for infinite case
k2=(~isfinite(N))&(p>0)&(rem(x,1)==0)&(x>=1);
if any(k2)
    %For each k1, create harmonic number vector
    %(Should be able to do it without a for-loop)
    kin=find(k2);
    for i=1:length(k2)
    hmvn(i)=sum(1./((1:floor(x(kin(i)))).^(1+p(kin(i)))));   
    end
    %Use zeta function
    y(k2)=hmvn./zeta(1+p(k2));
end


%Catch round-off
y(y>1)=1;

end