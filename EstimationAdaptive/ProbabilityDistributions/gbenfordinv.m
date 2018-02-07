function x = gbenfordinv(p,b,n)
%GBENFORDINV Inverse of the generalized Benford cdf
%   X = GBENFORDINV(P,B,N) returns the inverse cdf of the generalized 
%   Beford distribution of base B, at the N'th digit, at the values in P.
%
%   Since the Generalized Benford distribution is discrete, GBENFORDINV
%   returns the least integer X, such that the Generalized Benford cdf
%   evaluated at X, equals or exceeds P.
%
%   Type: Discrete, bounded, {0,...,B-1}
%   Restrictions:
%        N >= 2         (N integer)  [Digit Position]
%        B >= 2         (B integer)  [Base]
%
%   See also GBENFORDPDF, GBENFORDCDF, GBENFORDSTAT, GBENFORDFIT,
%            GBENFORDLIKE, GBENFORDRND, GBENFORDSF, GBENFORDHAZ

%   Mike Sheppard
%   Last Modified: 15-Dec-2011


if nargin < 3
    error('gbenfordinv:TooFewInputs',...
        'Requires at least three input argument.');
end


[errorcode p b n] = distchck(3,p,b,n);

if errorcode > 0
    error('gbenfordinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end



%Initialize X to 0.
if isa(p,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (b>=2 & b<Inf & b==round(b)) & (n>=2 & n<Inf & n==round(n));
okvar = (0 < p & p < 1);
k=(okparam & okvar);

x(~k)=NaN;
x(okparam & p==0)=1;
x(okparam & p==1)=b(okparam & p==1)-1; %one less than base

k=find(k);
cumdist = x;
cumdist(k) = gbenfordpdf(0,b(k),n(k));
count = 0;
k = k(cumdist(k) < p(k));

while ~isempty(k)
    x(k)=x(k)+1;
    count=count+1;
    cumdist(k)=cumdist(k)+gbenfordpdf(count,b(k),n(k));
    k=k(cumdist(k)<p(k));
end



end

