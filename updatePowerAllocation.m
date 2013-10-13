function An = updatePowerAllocation(K, Q, I, P, A, S, closures)
  An = A;
  for l = 1 : K
    for q = 1 : Q
      sub = S((l - 1) * Q + q, :);
      alloc = A((l - 1) * Q + q, :);
      direc = zeros(1, K * I);
      relaxH = 0;
      for k = closures(l, :)
      	if k == 0
          continue;
        end
        for i = 1 : I
          a = alloc((k - 1) * I + i);
          direc((k - 1) * I + i) = a - sub((k - 1) * I + i) * a / 23;
          if direc((k - 1) * I + i) > relaxH
          	relaxH = direc((k - 1) * I + i);
          end
        end
      end
      relaxL = 0;
      relax = (relaxL + relaxH) / 2;
      proj = waterFilling(K, I, closures(l, :), direc, relax, 1e-6);
      while abs(sum(proj) - P) / P > 1e-6
        if sum(proj) > P
          relaxL = relax;
        elseif sum(proj) < P
          relaxH = relax;
        end
        relax = (relaxL + relaxH) / 2;
        proj = waterFilling(K, I, closures(l, :), direc, relax, 1e-6);
      end
      An((l - 1) * Q + q, :) = proj;
    end
  end
  return

function proj = waterFilling(K, I, closure, direc, relax, epsilon)
  proj = zeros(1, K * I);
  for k = closure
  	if k == 0
      continue;
    end
  	for i = 1 : I
      j = (k - 1) * I + i;
      proj(j) = (direc(j) - relax <= epsilon) ? epsilon : direc(j) - relax;
  	end
  end
  return
