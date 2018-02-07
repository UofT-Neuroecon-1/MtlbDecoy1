function r = logrnd(mu,sigma,varargin)
%LOGRND Random arrays from the Logistic Distribution
%   R = LOGRND(MU,SIGMA) returns an array of random numbers chosen from the
%   Logistic Distribution with mean U and scale parameter SIGMA.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   The size of RND is the common size of the input arguments if both are
%   arrays.
%
%   R = LOGRND(MU,SIGMA,M,N,...) or R = LOGRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011



if nargin < 2
    error('logrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use loginv
r = loginv(rand(varargin{:}),mu,sigma);

end