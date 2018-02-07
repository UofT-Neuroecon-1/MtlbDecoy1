function r = yulernd(a,varargin)
%YULERND Random arrays from the Yule–Simon distribution
%   R = YULERND(A) returns an array of random numbers from the Yule-Simon
%   distribution with shape parameter A.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.

%   R = YULERND(A,M,N,...) or R = BBINORND(A,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 6-Apr-2011

if nargin < 1
    error('yulernd:TooFewInputs','Requires at least one input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Exponential-geometric mixture
w = exprnd(1./a,varargin{:});
r = geornd(exp(-w),varargin{:});


end

