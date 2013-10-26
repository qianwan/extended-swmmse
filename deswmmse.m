clear;
K = 1;
M = 4;
N = 2;
Q = 10;
I = 20;
SNRdB = 5;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusters = zeros(K, 1);
r = 2000;
if K == 1
    clusters = 0;
elseif K == 4
    clusters = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
end
closures = findClusterClosures(clusters, r * 0.9);
[bss, ues] = brownian(K, Q, I, clusters, r / sqrt(3));
H = generateMIMOChannel(K, Q, M, bss, I, N, ues, 2);
L = generateLambdas(K, Q, M, I, N, P, H, SNR, 2);
[V, A] = generateRandomTxVector(K, Q, M, I, P, closures);
[U, W, R, obj] = updateSWMmseVariables(K, Q, M, I, N, H, V);
reserve = 1e-6;
method = 'bcd';
numUnitIters = 3;

numAgents = 20;
LL = zeros(size(L, 1), size(L, 2) * numAgents);
LL(:, 1 : size(L, 2)) = L;
for i = 1 : size(LL, 1)
    for j = size(L, 2) + 1 : size(LL, 2)
        LL(i, j) = LL(i, j) + LL(i, j) * (rand - 0.5);
    end
end
CR = 0.5;
dw = 1;

maxIterations = 10;
numIterations = 0;
while numIterations < maxIterations
    numIterations = numIterations + 1;
    for i = 1 : numAgents
        otherAgents = setdiff(1 : numAgents, i);
        a = otherAgents(ceil(numel(otherAgents) * rand));
        otherAgents = setdiff(1 : numAgents, [i, a]);
        b = otherAgents(ceil(numel(otherAgents) * rand));
        otherAgents = setdiff(1 : numAgents, [i, a, b]);
        c = otherAgents(ceil(numel(otherAgents) * rand));
        Li = LL(:, (i - 1) * size(L, 2) + 1 : i * size(L, 2));
        [Ri, servBSsi, obji] = deSWMmseUnit(K, Q, M, I, N, H, U, V, W, Li, P, numUnitIters, reserve, method);
        La = LL(:, (a - 1) * size(L, 2) + 1 : a * size(L, 2));
        Lb = LL(:, (b - 1) * size(L, 2) + 1 : b * size(L, 2));
        Lc = LL(:, (c - 1) * size(L, 2) + 1 : c * size(L, 2));
        Rr = ceil(K * Q * rand);
        Rc = ceil(K * I * rand);
        Ln = Li;
        for q = 1 : K * Q
            for i = 1 : K * I
                ri = rand;
                if ri < CR || (q == Rr && i == Rc)
                    Ln(q, i) = La(q, i) + dw * (Lb(q, i) - Lc(q, i));
                else
                    Ln(q, i) = Li(q, i);
                end
            end
        end
        [Rn, servBSsn, objn] = deSWMmseUnit(K, Q, M, I, N, H, U, V, W, Ln, P, numUnitIters, reserve, method);
        if objn > obji
            LL(:, (i - 1) * size(L, 2) + 1 : i * size(L, 2)) = Ln;
        end
    end
end
