function C = unmixISTA(mes, endmemb, lambda)
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

    BTB = endmemb'*endmemb;
    gradL = @(C, BTs) BTB*C - BTs;
    T = @(x,t) (abs(x) - lambda * t) .* sign(x);

    % Optimize each spectrum in turn
    for i = 1:n
        S = mes(:,i);
        c = C(:,i);
        cLast = zeros(k,1);

        t = 1;
        BTs = endmemb'*S;
        while vecnorm(c-cLast,2) > 0.01
            cLast = c;
            step = 0.01/sqrt(t);
            
            c = T(c - 2 * step * gradL(c, BTs), step);
            c(c < 0) = 0;
            t = t + 1;
        end

        % Save the optimized abundance vector
        C(:,i) = c;
    end

end