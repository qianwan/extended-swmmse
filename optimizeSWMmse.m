function X = optimizeSWMmse(K, Q, M, I, J, D, V, L, P)
    X = V;
    for k = 1 : K
        lambda = L(k);
        while true
            % fprintf(2, '%d', k);
            T = zeros(Q, 1);
            offset = (k - 1) * Q * M;
            V(offset + 1 : offset + Q * M, :) = X(offset + 1 : offset + Q * M, :);
            for q = 1 : Q
                [C, A] = cvector(Q, M, I, D, J, X, lambda, k, q);
                for i = 1 : I
                    if A(i) == 0
                        rowOffset = (k - 1) * Q * M + (q - 1) * M;
                        colOffset = (k - 1) * I + i;
                        X(rowOffset + 1 : rowOffset + M, colOffset) = zeros(M, 1);
                    end
                end
                if nnz(A) == 0
                    T(q) = 0;
                    continue;
                end
                rowOffset = (k - 1) * Q * M + (q - 1) * M;
                colOffset = (q - 1) * M;
                Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                rho = spectralRadius(Jkq);
                miuLow = 0;
                miuHigh = upperBoundOfMiu(I, P, A, C);
                miu = (miuLow + miuHigh) / 2;
                delta = zeros(I, 1);
                rowOffset = (k - 1) * Q * M + (q - 1) * M;
                colOffset = (q - 1) * M;
                Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                while true
                    miu = (miuLow + miuHigh) / 2;
                    deltaLow = zeros(I, 1);
                    deltaHigh = upperBoundOfDelta(A, C, I, rho, lambda, miuHigh);
                    for i = 1 : I
                        if A(i) == 0
                            continue;
                        end
                        while true
                            delta(i) = (deltaLow(i) + deltaHigh(i)) / 2;
                            target = bisectionTarget(Jkq, C(:, i), delta(i), lambda, miu);
                            if target < 1
                                deltaLow(i) = delta(i);
                            else
                                deltaHigh(i) = delta(i);
                            end
                            if abs(deltaLow(i) - deltaHigh(i)) < 1e-6
                                break;
                            end
                        end
                    end
                    power = 0;
                    for i = 1 : I
                        if A(i) == 0
                            continue;
                        end
                        power = power + 1 / delta(i)^2;
                    end
                    if power < P
                        miuHigh = miu;
                    else
                        miuLow = miu;
                    end
                    if abs(miuLow - miuHigh) < 1e-6
                        break;
                    end
                end
                for i = 1 : I
                    if A(i) == 0
                        continue;
                    end
                    c = C(:, i);
                    v = (Jkq + (lambda * delta(i) / 2 + miu) * eye(M)) \ c;
                    rowOffset = (k - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    X(rowOffset + 1 : rowOffset + M, colOffset) = v;
                end
                T(q) = miu;
            end
            if checkSWMmseConverged(Q, M, I, T, P, X, D, J, k, lambda) == true
                break;
            end
        end
        % fprintf(2, '\n');
    end
    return

function y = checkSWMmseConverged(Q, M, I, T, P, V, D, J, k, lambda)
    y = true;
    for q = 1 : Q
        miu = T(q);
        [C, A] = cvector(Q, M, I, D, J, V, lambda, k, q);
        rowOffset = (k - 1) * Q * M + (q - 1) * M;
        colOffset = (q - 1) * M;
        Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
        for i = 1 : I
            rowOffset = (k - 1) * Q * M + (q - 1) * M;
            colOffset = (k - 1) * I + i;
            v = V(rowOffset + 1 : rowOffset + M, colOffset);
            if norm(v) == 0
                if norm(C(:, i)) <= lambda / 2
                    continue;
                else
                    y = false;
                    return;
                end
            else
                delta = 1 / norm(v);
                target = bisectionTarget(Jkq, C(:, i), delta, lambda, miu);
                if abs(target - 1) < 1e-6
                    continue;
                else
                    y = false;
                    return;
                end
            end
        end
        power = 0;
        for i = 1 : I
            rowOffset = (k - 1) * Q * M + (q - 1) * M;
            colOffset = (k - 1) * I + i;
            v = V(rowOffset + 1 : rowOffset + M, colOffset);
            power = power + dot(v, v);
        end
        if abs(miu * (P - power)) > 1e-6
            y = false;
            return
        end
    end
    return

function [C, A] = cvector(Q, M, I, D, J, V, lambda, k, q)
    C = zeros(M, I);
    for i = 1 : I
        rowOffset = (q - 1) * M;
        colOffset = (k - 1) * I + i;
        d = D(rowOffset + 1 : rowOffset + M, colOffset);
        c = d;
        for p = 1 : Q
            if p == q
                continue;
            end
            rowOffset = (k - 1) * Q * M + (q - 1) * M;
            colOffset = (p - 1) * M;
            Jkp = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
            rowOffset = (k - 1) * Q * M + (p - 1) * M;
            colOffset = (k - 1) * I + i;
            v = V(rowOffset + 1 : rowOffset + M, colOffset);
            c = c - Jkp * v;
        end
        C(:, i) = c;
    end
    A = zeros(I, 1);
    for i = 1 : I
        c = C(:, i);
        if norm(c, 2) > lambda / 2
            A(i) = i;
        end
    end
    return

function miuHigh = upperBoundOfMiu(I, P, A, C)
    maxc = 0;
    for i = 1 : I
        if A(i) == 0
            continue;
        end
        c = C(:, i);
        if norm(c) > maxc
            maxc = norm(c);
        end
    end
    miuHigh = (P / nnz(A))^(-0.5) * maxc * 1.1;
    return

function deltaHigh = upperBoundOfDelta(A, C, I, rho, lambda, miuHigh)
    deltaHigh = zeros(I, 1);
    for i = 1 : I
        if A(i) == 0
            continue;
        end
        deltaHigh(i) = (rho + miuHigh) / (norm(C(:, i)) - lambda / 2) * 1.1;
    end
    return

function rho = spectralRadius(A)
    rho = max(abs(eig(A)));
    return

function ret = bisectionTarget(Jkq, ciq, delta, lambda, miu)
    ret = delta * norm((Jkq + (lambda * delta / 2 + miu) * eye(size(Jkq))) \ ciq);
    return
