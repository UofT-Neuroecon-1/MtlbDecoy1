function r = gbenfordrnd(b,n,varargin)
%GBENFORDRND Random arrays from the generalized Benford distribution
%   R = GBENFORDRND(B,N) returns an array of random numbers chosen from the
%   Generalized Benford distribution with base B at the N'th digit
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = GBENFORDRND(B,N,MM,NN,...) or R = GBENFORD(B,N,[MM,NN,...]) 
%   returns an MM-by-NN-by-... array.
%
%   See also GBENFORDPDF, GBENFORDCDF, GBENFORDINV, GBENFORDSTAT, 
%            GBENFORDFIT, GBENFORDLIKE
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 2
    error('gbenfordrnd:TooFewInputs','Requires at least two input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Use gbenfordinv
r=gbenfordinv(rand(varargin{:}),b,n);

end