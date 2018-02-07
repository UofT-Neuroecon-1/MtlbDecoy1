function y = vanhcdf(x,a,pa,b,pb)
%VANHCDF Van Houtum cumulative distribution
%   Y = VANHCDF(A,PA,B,PB) returns the Van Houtum cumulative distribution with
%   parameters P(X=A)=PA, P(X=B)=PB and A<=B
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 3-Jun-2011


if nargin < 5
    error('vanhcdf:TooFewInputs',...
          'Requires five input arguments.'); 
end

[errorcode x a pa b pb] = distchck(5,x,a,pa,b,pb);

if errorcode > 0
    error('vanhcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(pa,'single') || ...
        isa(b,'single') || isa(pb,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


y(x<a)=0;
y(x==a)=pa(x==a);
y(x>=b)=1;


k=(x>a)&(x<b)&(~mod(x,1)&pa>=0&pa<=1&pb>=0&pb<=1);
if any(k)
    y(k)=pa(k)+floor(x(k)-a(k)).*((1-pa(k)-pb(k))./(b(k)-a(k)-1));
end

y(a>b) = NaN;

end