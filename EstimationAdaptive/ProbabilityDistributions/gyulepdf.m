function y = gyulepdf(x,a,b)
%GYULEPDF Generalized Yule–Simon probability density function
%   Y = GYULEPDF(X,A,B) returns the Generalized Yule-Simon probability 
%   density function with shape parameter A and B, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, semi-bounded, {1,...,Inf}
%   Restrictions:
%      A > 0
%      0 <= B < 1
%
%   NOTE: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1. The two-parameter generalization of the
%   original Yule distribution replaces the beta function with an
%   incomplete beta function. 
%
%   See also GYULECDF, GYULEINV, GYULESTAT, GYULEFIT,
%            GYULELIKE, GYULERND, GYULESF, GYULEHAZ
%

%   Mike Sheppard
%   Last Modified 15-May-2012


if nargin ~= 3
    error('gyulepdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('gyulepdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a) & (0<=b & b<1);
okvar = (1 <= x) & (x < Inf) & (x==round(x));
ok=(okparam & okvar);    
y(~okparam)=NaN;
y(okparam & ~okvar)=0;


if any(ok)
    a = a(ok); b = b(ok); x = x(ok);
    %Let Incomplete Beta catch errors on x,a, and b, yielding NaNs as outputs
    y(ok) = (a./(1-b.^a)).*(betainc(1-b,x,a+1).*beta(x,a+1));
end


%Catch round off
y(y<0)=0;


end