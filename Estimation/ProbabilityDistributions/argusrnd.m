function r = argusrnd(chi,c,varargin)
%ARGUSRND Random arrays from the ARGUS distribution
%   R = ARGUSRND(CHI,C) returns an array of random numbers chosen from the
%   ARGUS distribution with curvature parameter CHI, and cut-off
%   parameter C.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = ARGUSRND(CHI,C,M,N,...) or R = ARGUSRND(CHI,C,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, bounded
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSINV, ARGUSSTAT, ARGUSFIT, ARGUSLIKE
%

%   Mike Sheppard
%   Last Modified 3-Jul-2011


if nargin < 2
    error('argusrnd:TooFewInputs','Requires at least two input arguments.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<chi & chi<Inf);
chi(~okparam)=NaN; c(~okparam)=NaN;

%Use argusinv to catch any errors
r = argusinv(rand(varargin{:}),chi,c);


end
