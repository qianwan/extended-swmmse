K = 4;
M = 4;
N = 2;
Q = 20;
I = 40;
SNRdB = 25;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusterLocations = [0 + 0j, 0 + 2000j, 2000 * cos(pi / 6) + 2000 * sin(pi / 6) * 1j, -2000 * cos(pi / 6) + 2000 * sin(pi / 6) * 1j];
clusterClosures = findClusterClosures(clusterLocations, 1000);
[bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, 2000 / sqrt(3));

numCases = 100;
totalSumRate = 0;
totalNumIterations = 0;
maxIterations = 1e6;
epsilon = 1e-3;
for i = 1 : numCases
  numIterations = 0;
  prev = 0.0;
  H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations, 2);
  V = generateRandomTxVector(K, Q, M, I, P, clusterClosures);
  [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V);
  while abs(prev - sum(R)) > epsilon
    prev = sum(R);
    numIterations = numIterations + 1;
    if numIterations > maxIterations
      numIterations = numIterations - 1;
      break;
    end
    mmse = updateMmseMMatrix(K, Q, M, I, N, H, U, W);
    V = iterateWMMSE(K, Q, M, I, N, mmse, P, H, W, U);
    [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V);
  end
  fprintf(2, 'Case #%d: R = %f, # = %d\n', i, sum(R), numIterations);
  totalSumRate = totalSumRate + sum(R);
  totalNumIterations = totalNumIterations + numIterations;
end
fprintf(2, 'Avg sum rate: %f\n', totalSumRate / numCases);
fprintf(2, 'Avg number of iterations: %f\n', totalNumIterations / numCases);
