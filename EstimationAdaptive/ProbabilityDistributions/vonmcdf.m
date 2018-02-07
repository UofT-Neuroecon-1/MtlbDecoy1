function y = vonmcdf(x,u,k)
%VONMCDF Cumulative Von-Mises probability  function
%   Y = VONMCDF(X,U,K) returns the cumulative Von-Mises probability  function
%   with mean U and concentration K.
%
%   Sums the integral of the series
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 3
    error('vonmcdf:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode x u k] = distchck(3,x,u,k);

if errorcode > 0
    error('vonmcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


%Two calculations:
%[1]: At u-pi (lower limit for each input value of U)
%[2]: At x    (input values of X)
%CDF = [2] - [1].
y=vonmseriesint(x,u,k)-vonmseriesint(u-pi,u,k);



end


function approxint=vonmseriesint(x,u,k)
%Turn all into vectors
sz=size(x); x=x(:); u=u(:); k=k(:);
j=1:1e4;  j=j(:)'; %Sum the first 1e4 terms, turn into column vector
%Turn all inputs into matrices
%Rows: Individual inputs
%Columns: Values of j
jM=repmat(j,numel(x),1);
xM=repmat(x,1,numel(j));
uM=repmat(u,1,numel(j));
kM=repmat(k,1,numel(j));
%Compute two temporary matrices for multiplying
besM=besseli(j,k);
sinM=sin(jM.*(xM-uM))./jM;
%Sum the series across second dimension and compute approximation for the
%integral in the correct shape of the inputs
sumM=sum(besM.*sinM,2); %Column Vector
approxint=(1/(2*pi)).*((x+((2./besseli(0,k)).*sumM)));
approxint=reshape(approxint,sz);
end
