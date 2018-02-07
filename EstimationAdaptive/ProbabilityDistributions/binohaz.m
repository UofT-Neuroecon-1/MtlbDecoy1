function h = binohaz(x,n,p)
%BINOHAZ Binomial hazard function
%   H=BINOHAZ(X,N,P) returns the hazard function of the binomial 
%   distribution with parameters N and P at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also BINOPDF, BINOCDF, BINOINV, BINOSTAT, 
%            BINOFIT, BINOLIKE, BINORND, BINOSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 3
    error('binohaz:TooFewInputs',...
        'Requires three input argument.');
end


try
    yt = binopdf(x,n,p);
    st = binosf(x,n,p);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('binohaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

