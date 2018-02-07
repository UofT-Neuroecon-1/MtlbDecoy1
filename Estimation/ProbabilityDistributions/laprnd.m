function r = laprnd(u,b,varargin)
%LAPRND Random arrays from the Laplace Distribution
%   R = LAPRND(U,B) returns an array of random numbers chosen from the
%   Laplace Distribution with mean U and scale parameter B
%
%   Type: Continuous, unbounded
%   Restrictions:
%     B>0
%
%   The size of RND is the common size of the input arguments if both are
%   arrays.
%
%   R = LAPRND(U,B,M,N,...) or R = LAPRND(U,B,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011



if nargin < 2
    error('laprnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={[1] [1]}; end %Scalar


%Use lapinv
r = lapinv(rand(varargin{:}),u,b);

end