function r = loggamrnd(a,b,u,varargin)
%LOGGAMRND Random arrays from the Log-Gamma Distribution
%   R = LOGGAMRND(A,B,U) returns an array of random numbers chosen from the
%   Log-Gamma Distribution with shape parameters A and B and location
%   parameter U.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A,B>0
%     U>=0
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = LOGGAMRND(A,B,U,M,N,...) or R = LOGGAMRND(A,B,U,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011



if nargin < 3
    error('loggamrnd:TooFewInputs',...
          'Requires at least three input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use loggaminv
r = loggaminv(rand(varargin{:}),a,b,u);

end