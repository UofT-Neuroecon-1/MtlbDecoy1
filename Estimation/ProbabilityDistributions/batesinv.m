function x = batesinv(p,n,a,b)
%BATESINV Inverse of the Bates cumulative distribution function
%   X = BATESINV(P,N,A,B) returns the inverse cumulative distribution
%   function of the Bates Distribution of the mean of N independent
%   and identically distributed random variables uniformly distributed
%   continuously from A to B
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%      A < B
%      N > 1  (integer)
%
%   NOTE: Loss of accuracy occurs when N is greater than 25.
%
%   See also BATESPDF, BATESCDF, BATESSTAT, BATESFIT,
%            BATESLIKE, BATESRND, BATESSF, BATESHAZ
%

%   Mike Sheppard
%   Last Modified: 15-Dec-2011

if (nargin<2)
    error('batesinv:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin==2
    a=0; b=1;
elseif nargin==3
    b=1;
end


[errorcode p n a b] = distchck(4,p,n,a,b);

if errorcode > 0
    error('batesinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    x = zeros(size(p),'single');
else
    x = zeros(size(p));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>1 & n==round(n));
okvar = (0 < p & p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=a(okparam & p==0);
x(okparam & p==1)=b(okparam & p==1);

if any(ok)
    p=p(ok); n=n(ok); a=a(ok); b=b(ok);
    %Transformation of Irwin-Hall Distribution
    x(ok) = irwinhallinv(p,n,a,b)./n;
end


end