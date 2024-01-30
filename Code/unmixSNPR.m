function C = unmixSNPR(mes, endmemb, lambda)
    % lambda is regularization weight
    % mes is an m x n matrix of the measured spectra
    % endmemb is an m x k matrix of the endmember specra
    % C is a k x n matrix of the unmixed endmember abundances
    % At the end, endmemb * c = mes

    n = size(mes,2);
    k = size(endmemb,2);

    % Initial guess is nonegative linear least squares
    C = zeros(k,n);
    for i = 1:n
        [C(:,i), ~] = lsqnonneg(endmemb,mes(:,i));
    end

    % Objective function
    % L = @(C, s) sum(endmemb*C - s' * log(endmemb*C)) + lambda*sum(C);
    % Gradient of objective
    gradL = @(C, s) lambda*ones(k,1) + sum(endmemb .* repmat(1-s./(endmemb*C), [1 k]))';
    % Hessian of objective
    % hessL = @(C, s) hessian(endmemb, s, C);

    % Optimize each spectrum in turn
    for i = 1:n
        S = mes(:,i);
        c = C(:,i);
        cLast = zeros(k,1);

        b = zeros(k,1); u = 0.9;
        t = 1;
        while vecnorm(c-cLast,2) > 0.01
            cLast = c;
            step = 0.01/sqrt(t);
            
            b = u*b + gradL(c,S);
            c = c - step * b;
            c(c < 0) = 0;
            t = t + 1;
        end

        % Save the optimized abundance vector
        C(:,i) = c;
    end

end

function H = hessian(B, s, c)
    k = size(c,2);
    H = zeros(k);
    for i = 1:k
        for j = 1:k
            H(i,j) = sum(B(:,i).*B(:,j).*s./(B*c).^2);
        end
    end
end