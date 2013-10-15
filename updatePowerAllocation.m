function An = updatePowerAllocation(K, Q, I, P, A, S, closures, delta)
  An = A;
  for l = 1 : K
    for q = 1 : Q
      sub = -S((l - 1) * Q + q, :);
      alloc = A((l - 1) * Q + q, :);
      direc = zeros(1, K * I);
      relaxH = 0;
      for k = closures(l, :)
      	if k == 0
          continue;
        end
        for i = 1 : I
          a = alloc((k - 1) * I + i);
          direc((k - 1) * I + i) = a - sub((k - 1) * I + i) * a * delta;
          if direc((k - 1) * I + i) > relaxH
          	relaxH = direc((k - 1) * I + i);
          end
        end
      end
      relaxL = 0;
      relaxH = 100;
      relax = (relaxL + relaxH) / 2;
      epsilon = 1e-6;
      proj = waterFilling(K, I, closures(l, :), direc, relax, epsilon);
      if (sum(waterFilling(K, I, closures(l, :), direc, relaxH, epsilon)) > P) && (sum(waterFilling(K, I, closures(l, :), direc, relaxL, epsilon)) > P)
        fprintf(2, 'Stuck hereA\n');
      end
      while abs(sum(proj) - P) > 1e-5
        if sum(proj) > P
          relaxL = relax;
        elseif sum(proj) < P
          relaxH = relax;
        else
          break;
        end
        relax = (relaxL + relaxH) / 2;
        proj = waterFilling(K, I, closures(l, :), direc, relax, epsilon);
      end
      An((l - 1) * Q + q, :) = proj;
    end
  end
  return

function proj = waterFilling(K, I, closure, direc, relax, reverse)
  proj = zeros(1, K * I);
  for k = closure
  	if k == 0
      continue;
    end
  	for i = 1 : I
      j = (k - 1) * I + i;
      if direc(j) - relax <= reverse
        proj(j) = reverse;
      else
        proj(j) = direc(j) - relax;
      end
  	end
  end
  return
