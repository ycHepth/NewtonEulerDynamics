function S = crossMatrix(vector)
%CROSSMATRIX 
%   product cross product matrix to descript cross product in matrix form
%   THE INPUT vector must be a 3x1 vector
 z = vector(3,:);
 y = vector(2,:);
 x = vector(1,:);
 
 S = [ 0  -z  y;
       z   0  x;
      -y   x  0];
end

