function r = cauchyrnd(a,b,varargin)
%CAUCHYRND Random arrays from the Cauchy Distribution
%   R = CAUCHYRND(A,B) returns an array of random numbers chosen from the
%   Cauchy Distribution with location parameter A and scale parameter B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = CAUCHYRND(A,B,M,N,...) or R = CAUCHYRND(A,B,[M,N,...]) 
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, unbounded
%   Restrictions:
%        B > 0
%
%   See also CAUCHYPDF, CAUCHYCDF, CAUCHYINV, CAUCHYSTAT, 
%            CAUCHYFIT, CAUCHYLIKE
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 2
    error('cauchyrnd:TooFewInputs','Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Cauchy(0,1) is the ratio of two standard normals
X=normrnd(0,1,varargin{:})./normrnd(0,1,varargin{:});
%If X~Cauchy(0,1) then bX+a ~ Cauchy(a,b)
r=a+(b.*X);

end