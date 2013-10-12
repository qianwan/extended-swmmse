K = 3;
M = 2;
N = 2;
Q = 1;
I = 1;
SNRdB = 25;
SNR = 10^(SNRdB / 10);
P = SNR;
clusterLocations = [0 + 0j, 0 + 2000j, 0 - 2000j];
clusterClosures = [1; 2; 3];
[bsLocations, ueLocations] = brownian(2000 / sqrt(3), clusterLocations, Q, I);

numCases = 100;
totalSumRate = 0;
totalNumIterations = 0;
maxIterations = 1000;
epsilon = 1e-2;
for i = 1 : numCases
  numIterations = 0;
  prev = 0.0;
  H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations);
  V = generateRandomTxVector(K, Q, M, I, P, clusterClosures);
  [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V);
  while abs(prev - sum(W)) > epsilon
    prev = sum(W);
    numIterations = numIterations + 1;
    if numIterations > maxIterations
      break;
    end
    mmse = updateMmseMMatrix(K, Q, M, I, N, H, U, W);
    V = iterateWMMSE(K, Q, M, I, N, mmse, P, H, W, U);
    [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V);
  end
  fprintf('Case #%d: sum rate = %f, #iterations%d\n', i, sum(R), numIterations);
  totalSumRate = totalSumRate + sum(R);
  totalNumIterations = totalNumIterations + numIterations;
end
fprintf('Avg sum rate: %f\n', totalSumRate / numCases);
fprintf('Avg number of iterations: %f\n', totalNumIterations / numCases);
