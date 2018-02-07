function y = zipfpdf(x,p,N)
%ZIPFPDF Zipf probability density function
%   Y = ZIPFPDF returns the Zipf density function with
%   with parameter P and range N.
%
%   If N is not given default is Infinity, which is equivalent to the Zeta
%   Distribution
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 23-May-2011


if nargin < 2
    error('zipfpdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end

if nargin==2
    N=Inf;
end

[errorcode x p N] = distchck(3,x,p,N);

if errorcode > 0
    error('zipfpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(p,'single') || isa(N,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end



k1=isfinite(N)&(p>0)&(rem(x,1)==0)&(x<=N)&(x>=1);
if any(k1)
    %For each k1, create harmonic number vector
    %(Should be able to do it without a for-loop)
    kin=find(k1);
    for i=1:length(kin)
    hmv(i)=sum(1./((1:N(kin(i))).^(1+p(kin(i)))));
    end
    y(k1)=(x(k1).^(-1-p(k1)))./hmv;
end


k2=(~isfinite(N))&(p>0)&(rem(x,1)==0)&(x>=1);
if any(k2)
    %Use zeta function
    y(k2)=(x(k2).^(-1-p(k2)))./zeta(1+p(k2));
end

end