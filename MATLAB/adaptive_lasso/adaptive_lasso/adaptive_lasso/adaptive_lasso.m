function [beta_hat]=adaptive_lasso(x, y, true_beta, n_folds)
    [n, p] = size(x);
    
    % test the found beta using newly generated data
    new_x = randn(10 * n, p);
    new_y = new_x * true_beta;
    
    % stage one
    [mse, lambda] = cv_lasso(x, y, ones(p, 1), n_folds);
    beta_hat = my_lasso(x, y, lambda * ones(p, 1), 100);
    
    fprintf('test error at stage 1: %f\n', 0.5 * norm(new_y - new_x * beta_hat) ^2);
    
    % stage two
    nz = find(beta_hat ~= 0);
    new_penalty = 1./beta_hat(nz);
    
    [mse, lambda] = cv_lasso(x(:, nz), y, new_penalty, n_folds);
    beta_hat = my_lasso(x(:, nz), y, lambda * new_penalty, 100);
    
    fprintf('test error at stage 1: %f\n', 0.5 * norm(new_y - new_x(:, nz) * beta_hat) ^2);
    
end