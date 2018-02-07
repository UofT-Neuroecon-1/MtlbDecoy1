function [m,v] = bortanstat(alpha,n)
%BORTANSTAT Mean and variance for the Borel-Tanner Distribution
%   [M,V] = BORTANSTAT(ALPHA,N) returns the mean and variance of the
%   Borel-Tanner Distribution with shape parameters ALPHA and N.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Discrete, semi-bounded, {N,...,Inf}
%   Restrictions:
%        1 <= N         (N integer)
%        0 < ALPHA < 1
%
%   See also BORTANPDF, BORTANCDF, BORTANINV, BORTANFIT, 
%            BORTANLIKE, BORTANRND, BORTANSF, BORTANHAZ
%

%   Mike Sheppard
%   Last Modified: 17-Dec-2011


if nargin<2
    error('bortanstat:TooFewInputs',...
         'Requires at least two input argument.');
end


try
    % Return NaN if any arguments are outside of their respective limits.
    okparam = (0<alpha & alpha<1) & (n>=1 & n==round(n));

    m=n./(1-alpha);
    v=(n.*alpha)./((1-alpha).^3);
    
    m(~okparam)=NaN; v(~okparam)=NaN;
catch
    error('bortanpdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


end