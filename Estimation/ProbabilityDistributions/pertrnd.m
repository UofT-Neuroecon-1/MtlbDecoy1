function r = pertrnd(xmin,xmax,c,lambda,varargin)
%PERTRND Random arrays from the PERT Distribution
%   R = PERTRND(XMIN,XMAX,C,LAMBDA) returns an array of random numbers
%   chosen from the PERT Distribution with range XMIN to XMAX and maximum
%   at C with shape parameter LAMBDA.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%   R = PERTRND(XMIN,XMAX,C,LAMBDA,M,N,...) or 
%   R = PERTRND(XMIN,XMAX,C,LAMBDA,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 1-Jul-2011


if nargin < 4
    error('pertrnd:TooFewInputs',...
          'Requires at least four input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar


r = pertinv(rand(varargin{:}),xmin,xmax,c,lambda);


end
