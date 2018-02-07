function y = truncpdf(dist,x,params,xmin,max)
%TRUNCPDF Truncated Distribution
%   y = TRUNCPDF(DIST,X,XMIN,XMAX) returns the probability density function
%   of the truncated distribution of DIST, truncated from XMIN to XMAX, at
%   the values of X.

%   Mike Sheppard
%   Last Modified 13-May-2011

distpdf=[dist 'pdf'];
distcdf=[dist 'cdf'];

eval([distpdf '(' num2str(.3) ')'])



end
