function phat = nakafit(data)
%NAKALIKE Parameter estimate for the Nakagami Distribution
%   PHAT = NAKAFIT(DATA) returns the maximum likelihood estimate of the
%   parameters P of the Nakagami Distribution given DATA
%
%   PHAT(1)=U (Shape parameter)
%   PHAT(2)=W (Spread parameter)
%
%   For further options use MLE function directly

%   Mike Sheppard
%   Last Modified 19-Jun-2011

if nargin < 1
    error('nakafit:TooFewInputs',...
        'Requires one input argument.');
end

phat=mle(data(:),'distribution','nakagami');

end
