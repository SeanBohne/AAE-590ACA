function Tilda = cross_mat(vec)
% Determines the cross product 3x3 matrix given the input vector

Tilda = [0, -vec(3), vec(2); 
         vec(3), 0, -vec(1); 
         -vec(2), vec(1), 0];
end