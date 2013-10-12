function [BSs, UEs] = brownian(outerRadius, locations, Q, I)
  K = length(locations);
  BSs = zeros(K * Q, 1);
  UEs = zeros(I * I, 1);
  for i = 1 : K
    for j = 1 : Q
      x = 0;
      y = 0;
      while true
        x = (rand - 0.5) * 2 * outerRadius;
        y = (rand - 0.5) * 2 * outerRadius;
        mx = abs(x);
        my = abs(y);
        valid = true;
        if my > outerRadius * sin(pi / 3) || (sqrt(3) * mx + my > outerRadius * sqrt(3))
           valid = false;
        end
        if valid == true
          break;
        end
      end
      BSs((i - 1) * Q + j) = x + y * 1j + locations(i);
    end
    for j = 1 : I
      while true
        x = (rand - 0.5) * 2 * outerRadius;
        y = (rand - 0.5) * 2 * outerRadius;
        mx = abs(x);
        my = abs(y);
        valid = true;
        if my > outerRadius * sin(pi / 3) || (sqrt(3) * mx + my > outerRadius * sqrt(3))
           valid = false;
        end
        if valid == true
          break;
        end
      end
      UEs((i - 1) * I + j) = x + y * 1j + locations(i);
    end
  end
  return
