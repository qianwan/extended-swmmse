function [X, S] = optimizePSWMmseSubproblem(K, Q, M, I, N, A, closures, mmse, H, V, U, W, L)
  X = V;
  S = zeros(K * Q, K * I);
  while true
    C = updateCVectors(K, Q, M, I, N, closures, mmse, H, X, U, W);
    for k = 1 : K
      for i = 1 : I
        for l = closures(k, :)
          if l == 0
            continue;
          end
          for q = 1 : Q
            rowOffset = (l - 1) * Q * M + (q - 1) * M;
            colOffset = (k - 1) * I + i;
            c = C(rowOffset + 1 : rowOffset + M, colOffset);
            offset = (l - 1) * Q * M + (q - 1) * M;
            mm = mmse(offset + 1 : offset + M, offset + 1 : offset + M);
            a = A((l - 1) * Q + q, (k - 1) * I + i);
            lambda = L((l - 1) * Q + q, (k - 1) * I + i);
            [x, multiplier] = blockCoordinateDescent(c, mm, a, lambda);
            rowOffset = (l - 1) * Q * M + (q - 1) * M;
            colOffset = (k - 1) * I + i;
            X(rowOffset + 1 : rowOffset + M, colOffset) = x;
            S((l - 1) * Q + q, (k - 1) * I + i) = multiplier;
          end
        end
      end
    end
    if checkSubproblemConverged(K, Q, M, I, N, A, C, S, X, closures, mmse)
      break;
    end
  end
  return

function converged = checkSubproblemConverged(K, Q, M, I, N, A, C, L, S, X, closures, mmse)
  converged = true;
  for l = 1 : K
    for q = 1 : Q
      for k = closures(l, :)
        if k == 0
          continue;
        end
        for i = 1 : I
          rowOffset = (l - 1) * Q * M + (q - 1) * M;
          colOffset = (k - 1) * I + i;
          c = C(rowOffset + 1 : rowOffset + M, colOffset);
          lambda = L((l - 1) * Q + q, (k - 1) * I + i);
          rowOffset = (l - 1) * Q * M + (q - 1) * M;
          colOffset = (k - 1) * I + i;
          x = X(rowOffset + 1 : rowOffset + M, colOffset);
          if norm(c, 2) <= lambda / 2
            if norm(x, 2) ~= 0
              converged = false;
              return
            end
          else
            multiplier = S((l - 1) * Q + q, (k - 1) * I + i);
            offset = (l - 1) * Q * M + (q - 1) * M;
            mm = mmse(offset + 1 : offset + M, offset + 1 : offset + M);
            a = A((l - 1) * Q + q, (k - 1) * I + i);
            if abs(multiplier * (a - norm(x, 2)^2)) > 1e-6
              converged = false;
              return
            end
            sub = lambda * x / norm(x) + 2 * (multiplier * x + mm * x - c);
            if norm(sub, 2) > 1e-6
              converged = false;
              return
            end
          end
        end
      end
    end
  end
  return

function [x, multiplier] = blockCoordinateDescent(c, mm, a, lambda)
  if norm(c, 2) <= lambda / 2
    x = zeros(M, 1);
    multiplier = 0;
  else
    theta = 1 / sqrt(a);
    target = bisectionTarget(0, theta, lambda, mm, c);
    if target > 1
      miuLow = 0;
      miuHigh = norm(c) * theta;
      multiplier = (miuLow + miuHigh) / 2;
      theta = 1 / sqrt(a);
      target = bisectionTarget(multiplier, theta, lambda, mm, c);
      while abs(target - 1) > 1e-6
        if target > 1
          miuLow = multiplier;
        elseif target < 1
          miuHigh = multiplier;
        end
        multiplier = (miuLow + miuHigh) / 2;
        target = bisectionTarget(multiplier, theta, lambda, mm, c);
      end
    elseif target < 1
      tLow = 0;
      tHigh = (spectralRadius(mm) + norm(c, 2) * theta) / (norm(c, 2) - lambda / 2);
      multiplier = 0;
      theta = (tLow + tHigh) / 2;
      target = bisectionTarget(multiplier, theta, lambda, mm, c);
      while abs(target - 1) > 1e-6
        if target > 1
          tHigh = theta;
        elseif target < 1
          tLow = theta;
        end
        theta = (tLow + tHigh) / 2;
        target = bisectionTarget(multiplier, theta, lambda, mm, c);
      end
    end
    x = inv(mm + (lambda * theta / 2 + multiplier) * eye(size(mm))) * c;
  end
  return

function rho = spectralRadius(A)
  rho = max(abs(eig(A)));
  return

function r = bisectionTarget(multiplier, theta, lambda, mm, c)
  r = theta * norm(inv(mm + (lambda * theta / 2 + multiplier) * eye(size(mm))) * c)
  return

function C = updateCVectors(K, Q, M, I, N, closures, mmse, H, V, U, W)
  C = zeros(K * Q * M, K * I);
  for l = 1 : K
    for q = 1 : Q
      for k = 1 : K
        for i = 1 : I
          c = zeros(M, 1);
          for m = closures(l, : )
            for p = 1 : Q
              if m == l && p == q
                continue;
              end
              rowOffset = (l - 1) * Q * M + (q - 1) * M;
              colOffset = (m - 1) * Q * M + (p - 1) * M;
              mm = mmse(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
              offset = (m - 1) * Q * M + (p - 1) * M;
              v = V(offset + 1 : offset + M, (k - 1) * I + i);
              c = c - mm * v;
            end
          end
          w = W((k - 1) * I + i);
          rowOffset = (k - 1) * I * N + (i - 1) * N;
          colOffset = (l - 1) * Q * M + (q - 1) * M;
          h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + M);
          offset = (k - 1) * I * N + (i - 1) * N;
          u = U(offset + 1 : offset + N, : );
          c = c + w * h' * u;
          rowOffset = (l - 1) * Q * M + (q - 1) * M;
          colOffset = (k - 1) * I + i;
          C(rowOffset + 1 : rowOffset + M, colOffset) = c;
        end
      end
    end
  end
  return
