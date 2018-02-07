function p = ricecdf(x,s,sigma)
%RICECDF Rician cumulative distribution function
%   P = RICECDF(X,S,SIGMA) returns the cumulative distribution function of
%   the Rician Distribution with noncentrality parameter S and scale 
%   parameter SIGMA, at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, semi-bounded, (0, Inf)
%   Restrictions:
%      S >= 0
%      SIGMA > 0
%
%   See also RICEPDF, RICEINV, RICESTAT, RICEFIT,
%            RICELIKE, RICERND, RICESF, RICEHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin ~= 3
    error('ricecdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x s sigma] = distchck(3,x,s,sigma);

if errorcode > 0
    error('ricecdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(s,'single') || isa(sigma,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

% Return NaN for out of range parameters.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

%From addrice.m in statistics toolbox
%-----
x(x<0) = 0;
p = ncx2cdf((x./sigma).^2, 2, (s./sigma).^2);
%-----


%Catch round off
p(p<0)=0; p(p>1)=1;

end