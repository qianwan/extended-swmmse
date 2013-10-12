function H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations, model)
  H = zeros(K * I * N, K * Q * M);
  for i = 1 : size(H, 1)
    for j = 1 : size(H, 2)
      H(i, j) = (randn + randn * 1j) / sqrt(2);
    end
  end
  return
