function [m,v] = gasinstat(alpha,a,b)
%GASINSTAT Mean and variance for the generalized arcsine distribution.
%   [M,V] = GASINSTAT(ALPHA,A,B) returns the mean and variance of the
%   Generalized Arcsine Distribution with shape parameter ALPHA, on the
%   interval [A,B]. 
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Type: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        0 < ALPHA < 1
%
%   See also GASINPDF, GASINCDF, GASININV, GASINFIT, 
%            GASINLIKE, GASINRND, GASINSF, GASINHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if (nargin~=1)&&(nargin~=3)
    error('gasininv:TooFewInputs',...
          'Requires either one or three input arguments.'); 
end

if nargin == 1
    a = 0;
    b = 1;
end

[errorcode alpha a b] = distchck(3,alpha,a,b);

if errorcode > 0
    error('gasinstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

[m,v] = betastat(1-alpha,alpha);

%Transform
m=a+m.*(b-a);
v=v.*(b-a).^2;

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (0<alpha & alpha<1);
m(~okparam)=NaN;
v(~okparam)=NaN;

end