L = 4.0;

// À adapter en fonction du raffinement de maillage souhaité 
h = 0.1;

Point(1) = {-L, -L, 0, h};
Point(2) = {-L, L, 0, h};
Point(3) = {L, -L, 0, h};
Point(4) = {L, L, 0, h};
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};

Physical Line(1) = {1, 2, 3, 4};

Line Loop(1) = {1:4}; 
Plane Surface(1) = {1}; 

Physical Surface(1) = {1};

