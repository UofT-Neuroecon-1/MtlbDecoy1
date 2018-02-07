function [m,v] = kumarstat(a,b)
%KUMARSTAT Mean and variance for the Kumaraswamy Distribution
%   [M,V] = KUMARSTAT(A,B) returns the mean and variance for the
%   Kumaraswamy Distribution with shape parameters A and B.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, bounded
%   Restrictions:
%      A , B > 0
%
%   Note: The notation for the Kumaraswamy distribution used here is 
%   F(X<x) = 1 - (1-x^A)^B  
%
%   See also KUMARPDF, KUMARCDF, KUMARINV, KUMARFIT, KUMARLIKE, KUMARRND
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011


if nargin < 2
    error('kumarstat:TooFewInputs',...
          'Requires two input argument.');
end

% Return NaN for out of range parameters.
a(a<0)=NaN; b(b<0)=NaN;

try
    m1=b.*beta(1+(1./a),b);
    m2=b.*beta(1+(2./a),b);
    m=m1;
    v=m2-(m1.^2);
catch
    error('kumarstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end



end