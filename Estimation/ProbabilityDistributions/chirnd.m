function r = chirnd(v,varargin)
%CHIRND Random arrays from Chi Distribution.
%   R = CHIRND(V) returns an array of random numbers chosen from the
%   Chi Distribution with V degrees of freedom.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = CHIRND(V,M,N,...) or R = CHIRND(V,[M,N,...]) returns an
%   M-by-N-by-... array. 
%
%   Type: Continuous, semi-bounded
%
%   See also CHIPDF, CHICDF, CHIINV, CHISTAT, CHIFIT, CHILIKE
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011

if nargin < 1
    error('chirnd:TooFewInputs','Requires at least one input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Use chiinv
r = chiinv(rand(varargin{:}),v);

end

