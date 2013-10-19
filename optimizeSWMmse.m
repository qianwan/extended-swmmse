function X = optimizeSWMmse(K, Q, M, I, N, J, D, V, U, W, L, P)
    X = V;
    count = 0;
    maxCount = 100;
    while true
        count = count + 1;
        if count > maxCount
            break;
        end
        fprintf(2, '.');
        if mod(count, 50) == 0
            fprintf(2, '\n');
        end
        T = zeros(K, Q);
        for k = 1 : K
            lambda = L(k);
            for q = 1 : Q
                [C, A] = cvector(K, Q, M, I, N, D, J, X, lambda, k, q);
                for i = 1 : I
                    if A(i) == 0
                        rowOffset = (k - 1) * Q * M + (q - 1) * M;
                        colOffset = (k - 1) * I + i;
                        X(rowOffset + 1 : rowOffset + M, colOffset) = zeros(M, 1);
                    end
                end
                if nnz(A) == 0
                    T(k, q) = 0;
                    continue;
                end
                miuLow = 0;
                miuHigh = upperBoundOfMiu(Q, M, I, P, A, C, k, q);
                miu = (miuLow + miuHigh) / 2;
                deltaLow = zeros(I, 1);
                deltaHigh = zeros(I, 1); % upperBoundOfDelta(Q, M, A, C, J, I, k, q, lambda, miuHigh);
                delta = zeros(I, 1); % (deltaLow + deltaHigh) / 2;
                rowOffset = (k - 1) * Q * M + (q - 1) * M;
                colOffset = (q - 1) * M;
                Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                while true
                    miu = (miuLow + miuHigh) / 2;
                    deltaHigh = upperBoundOfDelta(Q, M, A, C, J, I, k, q, lambda, miuHigh);
                    for i = 1 : I
                        if A(i) == 0
                            continue;
                        end
                        c = C(:, i);
                        while true
                            delta(i) = (deltaLow(i) + deltaHigh(i)) / 2;
                            target = bisectionTarget(Jkq, c, delta(i), lambda, miu);
                            if target < 1
                                deltaLow(i) = delta(i);
                            elseif target > 1
                                deltaHigh(i) = delta(i);
                            end
                            if abs((deltaLow(i) - deltaHigh(i)) / deltaHigh) < 1e-3
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
                    elseif power > P
                        miuLow = miu;
                    end
                    if abs(miuLow - miuHigh) < 1e-3
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
                T(k, q) = miu;
            end
        end
        if checkSWMmseConverged(K, Q, M, I, T, P, X)
            break;
        end
    end
    return

function y = checkSWMmseConverged(K, Q, M, I, T, P, V)
    y = true;
    for k = 1 : K
        for q = 1 : Q
            miu = T(k, q);
            power = 0;
            for i = 1 : I
                rowOffset = (k - 1) * Q * M + (q - 1) * M;
                colOffset = (k - 1) * I + i;
                v = V(rowOffset + 1 : rowOffset + M, colOffset);
                power = power + dot(v, v);
            end
            if abs(miu * (P - power)) > 1e-4
                y = false;
                return
            end
        end
    end
    return

function [C, A] = cvector(K, Q, M, I, N, D, J, V, lambda, k, q)
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
    A = zeros(1, I);
    for i = 1 : I
        c = C(:, i);
        if norm(c, 2) > lambda / 2
            A(i) = i;
        end
    end
    return

function miuHigh = upperBoundOfMiu(Q, M, I, P, A, C, k, q)
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
    miuHigh = (P / nnz(A))^(-0.5) * maxc;
    return

function deltaHigh = upperBoundOfDelta(Q, M, A, C, J, I, k, q, lambda, miuHigh)
    deltaHigh = zeros(I, 1);
    rowOffset = (k - 1) * Q * M + (q - 1) * M;
    colOffset = (q - 1) * M;
    Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
    rho = spectralRadius(Jkq);
    for i = 1 : I
        if A(i) == 0
            continue;
        end
        deltaHigh(i) = (rho + miuHigh) / (norm(C(:, i)) - lambda / 2);
    end
    return

function rho = spectralRadius(A)
    rho = max(abs(eig(A)));
    return

function ret = bisectionTarget(Jkq, ciq, delta, lambda, miu)
    ret = delta * norm((Jkq + (lambda * delta / 2 + miu) * eye(size(Jkq))) \ ciq);
    return
