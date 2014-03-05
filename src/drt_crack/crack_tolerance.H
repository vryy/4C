#ifndef CRACK_TOLERANCE_H
#define CRACK_TOLERANCE_H


// This in in DEGREE, not in radians
// if angle between an edge and the exact propagation angle is less
// than this value, then the crack is assumed to propagate through this edge
#define MIN_PROP_ANGLE 0.1

// In order to find the correct node through which the crack propagates,
// we define a vector along the direction of propagation angle from crack tip
// The length of this vector should be larger than the appropriate dimension of elements
// so that this vector intersects one of the edges of the element
#define LENGTH_PROP_VECTOR 1.0

// If the difference between any two angle (in radians) are less than this value, then both angles
// are considered to be equal
#define ANGLE_TOL_ZERO 1e-5

#endif
