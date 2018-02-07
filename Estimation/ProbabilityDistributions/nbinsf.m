function s = nbinsf(x,r,p)
%NBINSF Negative binomial survival function.
%   S=NBINSF(X,R,P) returns the survival function of the negative binomial
%   distribution with parameters R and P at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also NBINPDF, NBINCDF, NBININV, NBINSTAT, NBINFIT, NBINLIKE, 
%            NBINRND, NBINHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 3, 
    error(message('nbinsf:TooFewInputs')); 
end 


try
    s = 1 - nbincdf(x,r,p);
catch
    error('nbinsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end