clear;
K = 4;
M = 4;
N = 2;
Q = 20;
I = 40;
SNRdB = 0;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusterLocations = zeros(1, K);
if K == 1
  clusterLocations = 0 + 0j;
elseif K == 4
  clusterLocations = [0 + 0j, ...
                      0 + 2000j, ...
                      2000 * cos(pi / 6) + 2000 * sin(pi / 6) * 1j, ...
                      -2000 * cos(pi / 6) + 2000 * sin(pi / 6) * 1j];
end
closures = findClusterClosures(clusterLocations, 2100);
L = ones(K * I, 1) * 0.3;
[bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, 2000 / sqrt(3));

numCases = 1;
totalSumRate = 0;
totalNumIterations = 0;
maxIterations = 10;
epsilon = 1e-1;
for i = 1 : numCases
  numIterations = 0;
  prev = 0;
  [bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, 2000 / sqrt(3));
  H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations, 2);
  [V, A] = generateRandomTxVector(K, Q, M, I, P, closures);
  [U, W, R] = updatePSWMmseVariables(K, Q, M, I, N, H, V);
  while abs(prev - sum(R)) > epsilon
    prev = sum(R);
    numIterations = numIterations + 1;
    if numIterations > maxIterations
      numIterations = numIterations - 1;
      break;
    end
    mmse = updatePSWMmseMatrix(K, Q, M, I, N, H, U, W);
    [V, S] = optimizePSWMmseSubproblem(K, Q, M, I, N, A, closures, mmse, H, V, U, W, L);
    fprintf(2, 'update varaibles\n');
    [U, W, R] = updatePSWMmseVariables(K, Q, M, I, N, H, V);
    fprintf(2, 'update power allocation\n');
    A = updatePowerAllocation(K, Q, I, P, A, S, closures, (sum(R) - prev) * 0.1);
    fprintf(2, '  Sum rate @#%d in case#%d: %f\n', numIterations, i , sum(R));
  end
  fprintf(2, '->Case #%d: R = %f # = %d\n', i, sum(R), numIterations);
  totalSumRate = totalSumRate + sum(R);
  totalNumIterations = totalNumIterations + numIterations;
end
fprintf(2, 'Avg sum rate: %f\n', totalSumRate / numCases);
fprintf(2, 'Avg number of iterations: %f\n', totalNumIterations / numCases);
