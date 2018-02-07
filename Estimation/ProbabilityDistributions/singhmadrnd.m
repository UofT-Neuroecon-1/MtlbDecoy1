function r = singhmadrnd(q,a,b,varargin)
%SINGHMADRND Random arrays from the Singh-Maddala distribution
%   R = SINGHMADRND(Q,A,B) returns an array of random numbers chosen
%   from the Singh-Maddala distribution with shape parameters Q and A and
%   scale parameter B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = SINGHMADRND(Q,A,B,M,N,...) or R = SINGHMADRND(Q,A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   The Singh-Maddala distribution is a special case of the Generalized
%   Beta Prime Distribution
%   SINGHMAPDF(Q,A,B) ~ GBETAPRPDF(1,Q,A,B)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 3
    error('singhmadrnd:TooFewInputs',...
          'Requires at least four input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Special case of gbetaprrnd
r = gbetaprrnd(1,q,a,b,varargin{:});

end