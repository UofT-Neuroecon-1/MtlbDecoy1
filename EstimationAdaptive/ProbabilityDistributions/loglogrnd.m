function r = loglogrnd(a,b,varargin)
%LOGLOGRND Random arrays from the Log-logistic Distribution
%   R = LOGLOGRND(A,B) returns an array of random numbers chosen from the
%   Log-logistic Distribution with scale parameter A and shape parameter B
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A>0
%     B>0
%
%   The size of RND is the common size of the input arguments if both are
%   arrays.
%
%   R = LOGLOGRND(A,B,M,N,...) or R = LOGLOGRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 23-Jun-2011



if nargin < 2
    error('loglogrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={[1] [1]}; end %Scalar


%Use logloginv
r = logloginv(rand(varargin{:}),a,b);

end
