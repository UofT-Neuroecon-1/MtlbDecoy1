function y = moyalpdf(x,u,s)
%MOYALPDF Moyal probability density function
%   Y = MOYALPDF(X,U,S) returns the probability density function of the
%   Moyal Distribution with location parameter U and scale parameter S.
%
%   Type: Continuous, unbounded
%   Restrictions:
%      S > 0
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011



if nargin < 3
    error('moyalpdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x u s] = distchck(3,x,u,s);

if errorcode > 0
    error('moyalpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(u,'single') || isa(s,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

k=(s>0);
if any(k)
    xk=x(k); uk=u(k); sk=s(k); 
    sq2pi=2.506628274631000502415765284811045253006986740609938316629;
    z=(xk-uk)./sk;
    y(k)=exp((-z/2)-(1/2).*exp(-z))./(sq2pi.*sk);
end

% Return NaN for out of range parameters.
y(s<0)=NaN;

%Round off
y(y<0)=0;

end