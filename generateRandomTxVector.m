function [V, A] = generateRandomTxVector(K, Q, M, I, P, closures)
    V = zeros(K * Q * M, K * I);
    A = zeros(K * Q, K * I);
    for l = 1 : K
        closure = closures(l, :);
        numUEs = nnz(closure) * I;
        power = P / numUEs;
        for q = 1 : Q
            for k = closure
                if k == 0
                    continue;
                end
                for i = 1 : I
                    v = randn(M, 1) + randn(M, 1) * 1j;
                    v = v / norm(v, 2) * sqrt(power);
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    V(rowOffset + 1 : rowOffset + M, colOffset) = v;
                    A((l - 1) * Q + q, (k - 1) * I + i) = power;
                end
            end
        end
    end
    return
