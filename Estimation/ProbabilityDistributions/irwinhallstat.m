function [m,v]=irwinhallstat(n,a,b)
%IRWINHALLSTAT Mean and variance for the Irwin-Hall Distribution.
%   [M,V] = IRWINHALLSTAT(N,A,B) returns the mean and variance of the 
%   Irwin-Hall distribution, of the sum of n independent random variables
%   uniformly distributed from A to B.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Type: Continuous, bounded, (A*N,B*N)
%   Restrictions:
%        A < B
%        N >= 1        (integer)
%
%   Note: The Irwin-Hall Distribution is also known as the
%   Uniform Sum Distribution. For statistics of the mean of
%   N uniformly distributed random variables, use BATESSTATS. 
%
%   See also IRWINHALLPDF, IRWINHALLCDF, IRWINHALLINV, IRWINHALLFIT, 
%            IRWINHALLLIKE, IRWINHALLRND, IRWINHALLSF, IRWINHALLHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin < 1
    error('irwinhallstat:TooFewInputs',...
          'Requires at least one input argument.'); 
end

if nargin == 1
    a=0; b=1;
elseif nargin==2
    b=1;
end


[errorcode n a b] = distchck(3,n,a,b);

if errorcode > 0
    error('irwinhallstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

if isa(n,'single') || isa(a,'single') || isa(b,'single')
   m = zeros(size(a),'single');
else
   m = zeros(size(a));
end
v = m;

% Return NaN if any arguments are outside of their respective limits.
ok = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>=1) & (n==round(n));
n(~ok) = NaN;

m=n.*(a+b)/2;
v=n.*((b-a).^2)/12;


end