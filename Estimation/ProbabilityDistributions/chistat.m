function [m,v]= chistat(nu)
%CHISTAT Mean and variance for the Chi Distribution.
%   [M,V] = CHISTAT(NU) returns the mean and variance of the 
%   Chi distribution with NU degrees of freedom.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        V > 0
%
%   See also CHIPDF, CHICDF, CHIINV, CHIFIT, CHILIKE, CHIRND, 
%            CHISF, CHIHAZ
%

%   Mike Sheppard
%   Last Modified 15-Mar-2011

if nargin < 1, 
    error('chistat:TooFewInputs','Requires one input argument.'); 
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<nu & nu<Inf);

m=sqrt(2).*exp(gammaln((nu+1)/2)-gammaln(nu/2));
v=nu-m.^2;

m(~okparam)=NaN; v(~okparam)=NaN;

end