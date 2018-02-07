function y = vanhpdf(x,a,pa,b,pb)
%VANHPDF Van Houtum Distribution probability density function
%   Y = VANHPDF(A,PA,B,PB) returns the probability density function
%   for the Van Houtum Distribution, defined as:
%      P(X=A)=PA
%      P(X=B)=PB
%   with all discrete values inbetween X=A and X=B are equally probable
%
%   Type: Discrete, bounded
%   Restrictions:
%        A <= X <= B  (integers)
%        (PA, PB) > 0
%        PA+PB <= 1
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 3-Jun-2011


if nargin < 5
    error('vanhpdf:TooFewInputs',...
          'Requires five input arguments.'); 
end

[errorcode x a pa b pb] = distchck(5,x,a,pa,b,pb);

if errorcode > 0
    error('vanhpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(pa,'single') || ...
        isa(b,'single') || isa(pb,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end



k=(x==round(x))&(pa>=0)&(pa<=1)&(pb>=0)&(pb<=1); %valid criteria
if any(k&(x>a)&(x<b))
    k2=k&(x>a)&(x<b);
    y(k2)=(1-pa(k2)-pb(k2))./(b(k2)-a(k2)-1);
end

if any(k&(x==a))
    ka=k&(x==a); y(ka)=pa(ka);
end

if any(k&(x==b))
    kb=k&(x==b); y(kb)=pb(kb);
end


end