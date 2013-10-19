function [U, W, R] = updateSWMmseVariables(K, Q, M, I, N, H, V)
  U = zeros(K * I * N, 1);
  W = zeros(K * I, 1);
  R = zeros(K * I, 1);
  for k = 1 : K
    for i = 1 : I
      C = zeros(N);
      for l = 1 : K
        for j = 1 : I
          rowOffset = (k - 1) * I * N + (i - 1) * N;
          colOffset = (l - 1) * Q * M;
          h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
          rowOffset = (l - 1) * Q * M;
          colOffset = (l - 1) * I + j;
          v = V(rowOffset + 1 : rowOffset + Q * M, colOffset);
          hv = h * v;
          C = C + hv * hv';
        end
      end
      C = C + eye(N);
      rowOffset = (k - 1) * I * N + (i - 1) * N;
      colOffset = (k - 1) * Q * M;
      h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
      rowOffset = (k - 1) * Q * M;
      colOffset = (k - 1) * I + i;
      v = V(rowOffset + 1 : rowOffset + Q * M, colOffset);
      localHv = h * v;
      u = C \ localHv;
      offset = (k - 1) * I * N + (i - 1) * N;
      U(offset + 1 : offset + N, :) = u;
      W((k - 1) * I + i) = 1 / (1 - real(dot(localHv, u)));
      localHvvH = localHv * localHv';
      L = C - localHvvH;
      R((k - 1) * I + i) = log2(real(det(eye(N) + localHvvH / L)));
    end
  end
  return
