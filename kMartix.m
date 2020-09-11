function K = kMartix(vector)
%KMARTIX
%   presented Iw in matrix form with inertia vector PI
 z = vector(3,:);
 y = vector(2,:);
 x = vector(1,:);
 
 K = [x y z 0 0 0;
      0 x 0 y z 0;
      0 0 x 0 y z;];
end

