function r = kumarrnd(a,b,varargin)
%KUMARRND Random arrays from the Kumaraswamy Distribution
%   R = KUMARRND(A,B) returns an array of random numbers chosen from the
%   Kumaraswamy Distribution with shape parameters A and B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   Type: Continuous, bounded
%   Restrictions:
%      A , B > 0
%
%   Note: The notation for the Kumaraswamy distribution used here is 
%   F(X<x) = 1 - (1-x^A)^B  
%
%   See also KUMARPDF, KUMARCDF, KUMARINV, KUMARSTAT, KUMARFIT, KUMARLIKE
%


%   Mike Sheppard
%   Last Modified 26-Jun-2011



if nargin < 2
    error('kumarrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use kumarinv
r = kumarinv(rand(varargin{:}),a,b);

end