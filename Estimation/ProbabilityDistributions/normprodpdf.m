function y = normprodpdf(x,s1,s2)
%NORMPRODPDF Normal Product Distribution probability density function
%   Y = NORMPRODPDF(X,S1,S2) returns the probability density function of
%   the Normal Product Distribution of the product of two Normal
%   Distributions with standard deviations S1 and S2
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      s1 , s2 >0
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 3
    error('normprodpdf:TooFewInputs',...
          'Requires at three input argument.');
end

[errorcode, x,s1,s2] = distchck(3,x,s1,s2);

if errorcode > 0
    error('normprodpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(s1,'single') || isa(s2,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end


k=(s1>0 & s2>0);
if any(k)
    prods=s1(k).*s2(k);
    K0 = besselk(0,abs(x(k))./prods);
    y(k)=K0./(pi*prods);
end

% Return NaN for out of range parameters.
y(s1<0 | s2<0)=NaN;

%Round off
y(y<0)=0;

end