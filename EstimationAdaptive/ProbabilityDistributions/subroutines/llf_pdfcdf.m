
%==========================================================================

function nll = llf_pdfcdf(params, uncensData, censData, uncensFreq, censFreq, pdfArgs, cdfArgs, checkVals, lb, ub)
% Given function handles to a PDF and a CDF, evaluate the negative
% log-likelihood for PARAMS given DATA.

pdfFun = pdfArgs{1};
pdfAddArgs = pdfArgs(2:end);
cdfFun = cdfArgs{1};
cdfAddArgs = cdfArgs(2:end);

% Bounds checking needed when function called by fminsearch.  When called
% by fmincon, the last two args will not be there.
if nargin > 8
    if any(params<=lb | ub<=params)
        nll = Inf;
        return
    end
end

% Log-likelihood = log(PDF(uncensored values) + log(1-CDF(censored values))
%
% First, evaluate the specified PDF of the uncensored data.
paramsCell = num2cell(params);
pdfVals = feval(pdfFun, uncensData, paramsCell{:}, pdfAddArgs{:});

% Make sure returned pdf values are valid.
if checkVals
    if any(~isfinite(pdfVals))
        error(message('stats:mle:NonfinitePdfVal'));
    elseif any(pdfVals <= 0)
        error(message('stats:mle:NonpositivePdfVal'));
    end
end

% Compute negative log-likelihood from uncensored values, using
% frequencies.
nll = -sum(uncensFreq.*log(pdfVals));

% If there is censoring, evaluate the specified CDF of the censored data.
if ~isempty(censData)
    cdfVals = feval(cdfFun, censData, paramsCell{:}, cdfAddArgs{:});

    % Make sure returned cdf values are valid.
    if checkVals
        if any(~isfinite(cdfVals))
            error(message('stats:mle:NonfiniteCdfVal'));
        elseif any(cdfVals < 0)
            error(message('stats:mle:NegativeCdfVal'));
        elseif any(cdfVals >= 1)
            error(message('stats:mle:GTOneCdfVal'));
        end
    end

    % Update negative log-likelihood with censored values, using
    % frequencies.
    nll = nll - sum(censFreq.*log(1-cdfVals));
end