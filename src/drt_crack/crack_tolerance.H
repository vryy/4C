#ifndef CRACK_TOLERANCE_H
#define CRACK_TOLERANCE_H


// This in in DEGREE, not in radians
// if angle between an edge and the exact propagation angle is less
// than this value, then the crack is assumed to propagate through this edge
//---
// Results are found to be sensitive for this quantity   (sudhakar 04/14)
// For one test case, when refining the mesh, it was necessary to modify this value
//----
#define MIN_PROP_ANGLE 0.5

// In order to find the correct node through which the crack propagates,
// we define a vector along the direction of propagation angle from crack tip
// The length of this vector should be larger than the appropriate dimension of elements
// so that this vector intersects one of the edges of the element
#define LENGTH_PROP_VECTOR 1.0

// If the difference between any two angle (in radians) are less than this value, then both angles
// are considered to be equal
#define ANGLE_TOL_ZERO 1e-5

// If ALE displacement of nodes are less than this value, then the nodes are not moved
#define ALE_LEN_TOL 1e-10

// If any value is less than this, then it is considered to be zero
#define APPROX_ZERO 1e-12

// If angle between diagonal node in a hex element and  propagation angle is less than this value, then
// crack propagates through this node; One HEX is split into two WEDGE
#define TOL_DIAG_NODE 0.175

#endif
