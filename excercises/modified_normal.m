% Modified Normal Matrix Generator
% Input: target_function (as a symbolic expression), variables (symbolic), regularization parameter
% Output: Modified normal matrix N_modified

function N_modified = generate_modified_normal_matrix(target_function, variables, lambda)

    % Input Check
    if nargin < 3
        error('Missing arguments: target_function, variables, and lambda are required.');
    end
    
    % Symbolic Variables
    n_vars = length(variables); % Number of unknowns
    A = sym(zeros(0, n_vars)); % Initialize empty design matrix
    
    % Step 1: Parse the Target Function
    fprintf('Target Function: %s\n', char(target_function));

    % Extract Residual Terms and Smoothness Terms
    residual_terms = sym(0);
    regularization_terms = sym(0);

    % Decompose the target function into parts
    terms = children(target_function);
    for i = 1:length(terms)
        term = terms(i);
        if has(term, '^2') % Square terms represent residuals or smoothness
            if has(term, diff(variables)) % Smoothness term (differences between variables)
                regularization_terms = regularization_terms + term;
            else % Residual term
                residual_terms = residual_terms + term;
            end
        end
    end

    % Step 2: Generate Design Matrix from Residual Terms
    fprintf('Residual Terms: %s\n', char(residual_terms));
    residual_gradient = gradient(residual_terms, variables);
    A = jacobian(residual_gradient, variables);
    N_residual = simplify(A' * A); % Normal matrix from residuals

    % Step 3: Generate Smoothness Matrix from Regularization Terms
    fprintf('Regularization Terms: %s\n', char(regularization_terms));
    regularization_gradient = gradient(regularization_terms, variables);
    D = jacobian(regularization_gradient, variables);
    N_regularization = lambda * simplify(D' * D);

    % Step 4: Combine to Get the Modified Normal Matrix
    N_modified = simplify(N_residual + N_regularization);

    % Display Results
    fprintf('Residual Normal Matrix (N_residual):\n');
    disp(N_residual);
    
    fprintf('Regularization Normal Matrix (N_regularization):\n');
    disp(N_regularization);
    
    fprintf('Modified Normal Matrix (N_modified):\n');
    disp(N_modified);
end

% Example Usage
% Define symbolic variables
syms yA yB yC lambda real

% Input: Target Function
target_function = (yB - yA)^2 + (yC - yB)^2; % Smoothness penalty
variables = [yA, yB, yC];

% Call the function
N_modified = generate_modified_normal_matrix(target_function, variables, lambda);
