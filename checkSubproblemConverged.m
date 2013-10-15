
function converged = checkSubproblemConverged(K, Q, M, I, A, C, L, S, X, closures, mmse)
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
          lambda = L((k - 1) * I + i);
          rowOffset = (l - 1) * Q * M + (q - 1) * M;
          colOffset = (k - 1) * I + i;
          x = X(rowOffset + 1 : rowOffset + M, colOffset);
          if norm(c, 2) <= lambda / 2
            if norm(x, 2) ~= 0
              converged = false;
              return
            end
          elseif norm(x) > 0
            multiplier = S((l - 1) * Q + q, (k - 1) * I + i);
            offset = (l - 1) * Q * M + (q - 1) * M;
            mm = mmse(offset + 1 : offset + M, offset + 1 : offset + M);
            a = A((l - 1) * Q + q, (k - 1) * I + i);
            if abs(multiplier * (a - norm(x, 2)^2)) > 1e-6
              converged = false;
              return
            end
            sub = lambda * x / norm(x) + 2 * (multiplier * x + mm * x - c);
            if norm(sub, 2) > 1e-5
              converged = false;
              return
            end
          else
            converged = false;
            return
          end
        end
      end
    end
  end
  return
