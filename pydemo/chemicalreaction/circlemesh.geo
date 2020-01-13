// Gmsh project created on Thu Nov 28 14:40:07 2019
//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 75, 75};
//// For getting quad elements
//Recombine Surface {3, 1, 2};
