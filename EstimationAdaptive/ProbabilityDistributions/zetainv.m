function x = zetainv(y,p)
%ZETAINV Inverse of the Zeta cumulative distribution function
%   X = ZETAINV(Y,P) returns the inverse of the Zeta cumulative
%   distribution with parameter P
%

%   Mike Sheppard
%   Last Modified 30-May-2011



if nargin < 2
    error('zetainv:TooFewInputs',...
          'Requires at least two input arguments.'); 
end


try
    %Expand size if necessary
    y=y+zeros(size(p));
    p=p+zeros(size(y));
catch
    error('zetainv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


k=1:numel(p);
k1=find(p<=0);


%Initialize X to 0.
if isa(y,'single') || isa(p,'single')
    x=zeros(size(y),'single');
else
    x=zeros(size(y));
end

cumdist=x;
x(k1)=NaN;
k([k1(:)])=[];
if isempty(k), return; end

cumdist(k)=zetapdf(1,p(k));

count=1;

k=k(cumdist(k)<y(k));
while ~isempty(k)
    x(k)=x(k)+1;
    count=count+1;
    cumdist(k)=cumdist(k)+zetapdf(count,p(k));
    k=k(cumdist(k)<y(k));
end

%One off
x=x+1;

end