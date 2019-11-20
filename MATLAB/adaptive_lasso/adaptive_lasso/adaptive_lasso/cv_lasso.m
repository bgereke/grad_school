function [mse, lambda]=cv_lasso(x, y, penalty, n_folds)

    [n, p] = size(x);
    lambdas = [0.0001, 0.001, 0.01, 0.1, 1, 10];
    n_lambda = length(lambdas);
    
    idx = randperm(n);
    
    fold_size = n/n_folds;
    
    mse = zeros(n_lambda, 1);

    for l=1:n_lambda
        
        for fold=1:n_folds
            if fold ~= n_folds
                test_idx = idx((fold-1)* fold_size+1: fold*fold_size);
                train_idx = [idx(1: (fold-1)* fold_size), idx(fold*fold_size+1: end)];
            else
                test_idx = idx((fold-1)* fold_size+1: end);
                train_idx = idx(1: (fold-1)* fold_size);
            end
            
            test_x = x(test_idx, :);
            test_y = y(test_idx);
            
            train_x = x(train_idx, :);
            train_y = y(train_idx);
            
            cur_lambda = lambdas(l);
            
            beta_hat = my_lasso(train_x, train_y, cur_lambda * penalty,100);
            
            % prediction error on the test fold
            mse(l)= mse(l) + 0.5 * norm(test_x * beta_hat - test_y)^2;
        end
    end
    
    fprintf('lambda = %f\t mse=%f\n', [lambdas; mse']);
    
    [mse, idx]=min(mse);
    lambda=lambdas(idx);
end