function r = levyrnd(mu,sigma,varargin)
%LEVYRND Random arrays from the Lévy Distribution
%   R = LEVYRND(MU,SIGMA) returns an array of random numbers chosen from
%   the Lévy Distribution with location parameter MU and dispersion
%   parameter SIGMA. 
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = LEVYRND(MU,SIGMA,M,N,...) or R = LEVYRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%     SIGMA > 0
%
%   NOTE: The Lévy Distribution is also known as a van der Waals profile,
%   and if MU=0 it is a special case of the Inverse-Gamma Distribution.
%
%   See also LEVYPDF, LEVYCDF, LEVYINV, LEVYSTAT, 
%            LEVYFIT, LEVYLIKE, LEVYSF, LEVYHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012



if nargin < 2
    error('levyrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use levyinv
r = levyinv(rand(varargin{:}),mu,sigma);

end