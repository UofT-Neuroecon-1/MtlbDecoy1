function r = loglaprnd(u,b,varargin)
%LOGLAPRND Random arrays from the Log-Laplace Distribution
%   R = LOGLAPRND(U,B) returns an array of random numbers chosen from the
%   Log-Laplce Distribution with distribution parameters U and B. 
%   U and B are the mean and scale parameter, respectively, of the 
%   associated Laplace distribution.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.

%   Type: Continuous, unbounded
%   Restrictions:
%     B>0
%
%   The size of RND is the common size of the input arguments if both are
%   arrays.
%
%   R = LOGLAPRND(U,B,M,N,...) or R = LOGLAPRND(U,B,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011



if nargin < 2
    error('loglaprnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use loglapinv
r = loglapinv(rand(varargin{:}),u,b);

end