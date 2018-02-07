function x = bortaninv(p,alpha,n)
%BORTANINV Inverse of the Borel-Tanner cumulative distribution function
%   X = BORTANINV(P,ALPHA,N) returns the inverse cumulative distribution
%   function of the Borel-Tanner distribution with shape parameters ALPHA
%   and N, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Since the Borel-Tanner distribution is discrete, BORTANINV
%   returns the least integer X such that the Borel-Tanner cdf
%   evaluated at X, equals or exceeds P.
%
%   Distribution: Discrete, semi-bounded, {N,...,Inf}
%   Restrictions:
%        1 <= N         (N integer)
%        0 < ALPHA < 1
%
%   See also BORTANPDF, BORTANCDF, BORTANSTAT, BORTANFIT,
%            BORTANLIKE, BORTANRND, BORTANSF, BORTANHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011



if nargin<3
    error('bortancdf:TooFewInputs','Requires three input argument.');
end

[errorcode p alpha n] = distchck(3,p,alpha,n);

if errorcode > 0
    error('bortancdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(p,'single') || isa(alpha,'single') || isa(n,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<alpha & alpha<1) & (n>=1 & n==round(n));
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=n(okparam & p==0);
x(okparam & p==1)=Inf;

if any(ok)
    k=find(ok);
    cumdist=x;
    if isempty(k), return; end
    cumdist(k)=bortanpdf(0,alpha(k),n(k));
    count=0;
    k=k(cumdist(k)<p(k));
    while ~isempty(k)
        x(k)=x(k)+1;
        count=count+1;
        cumdist(k)=cumdist(k)+bortanpdf(count,alpha(k),n(k));
        k=k(cumdist(k)<p(k));
    end
end


end