function r = benfordrnd(b,varargin)
%BENFORDRND Random arrays from the Benford Distribution
%   R = BENFORDRND(B) returns an array of random numbers chosen from the
%   Benford distribution with base B. 
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BENFORDRND(B,M,N,...) or R = BENFORDRND(B,[M,N,...]) 
%   returns an M-by-N-by-... array.
%
%   See also BENFORDPDF, BENFORDCDF, BENFORDINV, BENFORDSTAT, 
%            BENFORDFIT, BENFORDLIKE
%

%   Mike Sheppard
%   Last Modified 13-Mar-2011


if nargin < 1
    error('benfordrnd:TooFewInputs','Requires at least one input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Could use IRWINHALLRND instead, but short enough to make it explicit
r=benfordinv(rand(varargin{:}),b);

end