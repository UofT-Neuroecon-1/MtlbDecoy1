function r = elrnd(p,b,varargin)
%ELRND Random arrays from the Exponential-logarithmic distribution
%   R = ELRND(P,B) returns an array of random numbers chosen from the
%   Exponential-logarithmic distribution with paramaters P (0<P<1)
%   and B (B>0)
%   The size of R is the common size of the inputs if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = ELRND(P,B,M,N,...) or R = ELRND(P,B,[M,N,...])
%   returns an M-by-N-by-... array.
%

%   Mike Sheppard
%   Last Modified 10-May-2011


if nargin < 2
    error('elrnd:TooFewInputs',...
          'Requires at least three input argument.'); 
end

if isempty(varargin), varargin={[1] [1]}; end %Scalar


sz=cell2mat(varargin);

%Using expgaminv
U=rand(sz);
r=(1/b).*log((1-p)./(1-(p.^U)));

end