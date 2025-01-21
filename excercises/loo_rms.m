% Leave-One-Out Cross Validation and RMSE Function with Weights
% Input: Dataset with t (time), y0 (observations), weight function
% Output: LOO predictions, RMSE, and weights

function [y_pred_loo, rmse, loo_weights] = loo_rms(t, y0, weight_function)

    % Input Check
    if nargin < 3
        error('Missing arguments: t, y0, weight_function are required.');
    end
    
    % Initialize Outputs
    n = length(t); % Number of observations
    y_pred_loo = zeros(n, 1); % LOO predictions
    errors = zeros(n, 1); % Errors for RMSE calculation
    loo_weights = cell(n, 1); % Store weights for each LOO iteration

    % Leave-One-Out Cross Validation
    for i = 1:n
        % Leave out observation i
        t_loo = t([1:i-1, i+1:end]);
        y0_loo = y0([1:i-1, i+1:end]);
        
        % Calculate weights for remaining observations
        weights = weight_function(abs(t_loo - t(i)));
        loo_weights{i} = weights; % Store the weights
        
        valid_weights = weights > 0;
        if sum(valid_weights) > 0
            y_pred_loo(i) = sum(weights(valid_weights) .* y0_loo(valid_weights)) / sum(weights(valid_weights));
        else
            y_pred_loo(i) = NaN; % No valid prediction
        end

        % Calculate error for LOO
        errors(i) = y0(i) - y_pred_loo(i);
    end

    % Root Mean Square Error (RMSE)
    rmse = sqrt(mean(errors(~isnan(errors)).^2));
end
