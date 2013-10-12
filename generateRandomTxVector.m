function V = generateRandomTxVector(K, Q, M, I, P, clusterClosures)
  V = zeros(K * Q * M, K * I);
  for m = 1 : K
    closures = clusterClosures(m, :);
    numUEs = nnz(closures) * I;
    power = P / numUEs;
    for q = 1 : Q
      for c = closures
        if c == 0
          continue;
        end
        for i = 1 : I
          v = randn(M, 1) + randn(M, 1) * 1j;
          v = v / norm(v, 2) * sqrt(power);
          rowOffset = (m - 1) * Q * M + (q - 1) * M;
          colOffset = (c - 1) * I + i;
          V(rowOffset + 1 : rowOffset + M, colOffset) = v;
        end
      end
    end
  end
  return
