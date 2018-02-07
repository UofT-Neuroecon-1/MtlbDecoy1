function y = zetacdf(x,p)
%ZETACDF Cumulative Zeta probability density function
%   Y = ZETACDF returns the cumulative Zeta density function with
%   with parameter P.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 25-May-2011


if nargin < 2
    error('zetacdf:TooFewInputs',...
        'Requires at least two input arguments.');
end


try
    %Expand size if necessary
    x=x+zeros(size(p));
    p=p+zeros(size(x));
catch
    error('zetacdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(p,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end


%Valid cases for infinite case
k2=(p>0)&(rem(x,1)==0)&(x>=1);
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


end