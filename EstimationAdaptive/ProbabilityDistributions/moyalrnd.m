function r = moyalrnd(u,s,varargin)
%MOYALRND Random arrays from the Moyal Distribution
%   R =MOYALRND(U,S) returns an array of random numbers chosen from the
%   Moyal Distribution  with location parameter U and scale parameter S.
%
%   Type: Continuous, unbounded
%   Restrictions:
%      U any real number
%      S>0
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = MOYALRND(U,S,M,N,...) or R = NAKARND(U,S,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 20-Jun-2011



if nargin < 2
    error('moyalstat:TooFewInputs',...
          'Requires two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use moyalinv
r = moyalinv(rand(varargin{:}),u,s);

end

