function [parmhat, parmci] = pascalfit(x, alpha,options)
%PASCALFIT Parameter estimates for Pascal Distribution data
%   PASCALFIT(X) returns the maximum likelihood estimates of the parameters
%   of the Pascal Distribution given the data in the vector, X.
%
%   [PARMHAT, PARMCI] = PASCAL(X,ALPHA) returns MLEs and 100(1-ALPHA)
%   percent confidence intervals given the data.  By default, the
%   optional parameter ALPHA = 0.05 corresponding to 95% conf. intervals.
%
%   [ ... ] = NBINFIT( ..., OPTIONS) specifies control parameters for the
%   numerical optimization used to compute ML estimates.  This argument can
%   be created by a call to STATSET.  See STATSET('nbinfit') for parameter
%   names and default values. (Note: PASCALFIT uses NBINFIT in its
%   computations)
%

%   Mike Sheppard
%   Last Modified 18-Jun-2011


% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
if nargin < 3 || isempty(options)
    options = statset('nbinfit');
else
    options = statset(statset('nbinfit'),options);
end
if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
p_int = [alpha/2; 1-alpha/2];

if min(size(x)) > 1
    error('stats:nbinfit:InvalidData','The first argument must be a vector.');
end

%Assume the minimum value N, then reiterate MLE of nbinfit until
%convergence
mx=min(x);
n=mx; ind=1;
while ind
[parmhat, parmci] = nbinfit(x-n, alpha, options);
q=round(parmhat(1));
if (q==n) | (q>mx), ind=0; end
n=q;
end

end
