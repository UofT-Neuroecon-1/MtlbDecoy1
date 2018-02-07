function x = zipfinv(y,p,N)
%ZIPFINV Inverse of the Zipf cumulative distribution function
%   X = ZIPFINV(Y,P,N) returns the inverse of the Zipf cumulative
%   distribution with parameter P and range N.
%
%   If N is not given default is Infinity, which is equivalent to the Zeta
%   Distribution
%
%   Mike Sheppard
%   Last Modified 30-May-2011



if nargin < 2
    error('zipfinv:TooFewInputs',...
          'Requires at least two input arguments.'); 
end

if nargin==2
    N=Inf;
end

[errorcode y p N] = distchck(3,y,p,N);

if errorcode > 0
    error('zipfinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end



k=1:numel(p);
k1=find(p<=0 | N~=round(N));


%Initialize X to 0.
if isa(y,'single') || isa(p,'single') || isa(N,'single')
    x=zeros(size(y),'single');
else
    x=zeros(size(y));
end

cumdist=x;
x(k1)=NaN;
k([k1(:)])=[];
if isempty(k), return; end

cumdist(k)=zipfpdf(1,p(k),N(k));

count=1;

k=k(cumdist(k)<y(k));
while ~isempty(k)
    x(k)=x(k)+1;
    count=count+1;
    cumdist(k)=cumdist(k)+zipfpdf(count,p(k),N(k));
    k=k(cumdist(k)<y(k));
end

%One off
x=x+1;

end
