// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 4
//
//  Built-in functions, holes in surfaces, annotations, entity colors
//
// -----------------------------------------------------------------------------

// As usual, we start by defining some variables:
lc = 0.3;
r = 2;
// We can use all the usual mathematical functions (note the capitalized first
// letters), plus some useful functions like Hypot(a, b) := Sqrt(a^2 + b^2):


// Then we define some points and some lines using these variables:

// Exterior corners of rectangle
Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, -8, 0, lc};
Point(4) = {0, -8, 0, lc};

// points to define circle
Point(5) = {5, -4, 0, lc}; // center
Point(6) = {5 - r, -4, 0, lc}; // right
Point(7) = {5 + r, -4, 0, lc}; // left

// exterior rectangle lines
Line(1) = {1, 2}; 
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// exterior rectangle curve loop
Curve Loop(5) = {1, 2, 3, 4};

// circle curves need 2 as arc must be less than pi
Circle(6) = {6, 5, 7};
Circle(7) = {7, 5, 6};
Curve Loop(8) = {6, 7};

// plane surface for interior circle
Plane Surface(1) = {8};

// plane surface for exterior rectangle with interior circle hole
Plane Surface(2) = {5, 8};

Physical Curve("HorEdges", 9) = {1, 3};
Physical Curve("VerEdges", 10) = {2, 4};
Physical Curve("Circle", 11) = {7, 6};

Physical Surface("PunchedDom", 3) = {2};
Physical Surface("Disc", 4) = {1};