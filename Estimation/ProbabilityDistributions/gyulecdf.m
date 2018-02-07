function p = gyulecdf(x,a,b)
%GYULECDF Generalized Yule–Simon cumulative distribution function
%   P = GYULECDF(X,A,B) returns the Generalized Yule-Simon cumulative 
%   distribution function with shape parameter A and B, at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
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
%   See also GYULEPDF, GYULEINV, GYULESTAT, GYULEFIT,
%            GYULELIKE, GYULERND, GYULESF, GYULEHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


% Superstars without talent?
% The Yule distribution controversy
% Laura Spierdijka and Mark Voorneveldb, c,1

if nargin ~= 3
    error('gyulecdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('gyulecdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize p to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a) & (0<=b & b<1);
okvar = (1 <= x) & (x < Inf) & (x==round(x));
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<=1)=0;
p(okparam & x==Inf)=1;


if any(ok)
    a = a(ok); b = b(ok); x = x(ok);
    %Let Incomplete Beta catch errors on x,a, and b, yielding NaNs as outputs
    p(ok) = 1-((a./(1-b.^a)).*(betainc(1-b,x+1,a).*beta(x+1,a)));
end


%Catch round off
p(p<0)=0; p(p>1)=1;


end