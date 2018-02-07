function r = baldnichrnd(F,p,varargin)
%BALDNICHRND Random arrays from the Balding-Nichols distribution
%   R = BALDNICHRND(F,P) returns an array of random numbers chosen from the
%   Balding-Nichols distribution with parameters F and p.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BALDNICHRND(F,P,M,N,...) or R = BALDNICHRND(F,P,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   The distribution assumes background allele frequency P the allele
%   frequencies, in sub-populations separated by Wright's F_ST F.
%
%   Type: Continuous, bounded
%   Restrictions:
%        0 < F, P < 1
%
%   Note: The Balding-Nichols distribution is a reparametrization of the
%   Beta distribution
%
%   See also BALDNICHPDF, BALDNICHCDF, BALDNICHINV, BALDNICHSTAT, 
%            BALDNICHFIT, BALDNICHLIKE
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin < 2
    error('baldnichrnd:TooFewInputs','Requires at least two input arguments.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<F & F<1) & (0<p & p<1);
F(~okparam)=NaN; p(~okparam)=NaN;

%Reparametrization of Beta Distribution
r = betarnd((1-F).*p./F,(1-F).*(1-p)./F,varargin{:});

end
