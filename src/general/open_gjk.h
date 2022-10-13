// Copyright 2022 Mattia Montanari, University of Oxford                          //
//                                                                                //
// This program is free software: you can redistribute it and/or modify it under  //
// the terms of the GNU General Public License as published by the Free Software  //
// Foundation, either version 3 of the License.

/* used for the collision detection */

#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

/**
 * @brief Macro that implements the CompareSign function (see paper).
 */
#define SAMESIGN( a, b ) ( (a>0) == (b>0) )


/**
 * @brief Structure of a body.
 */

struct bd {
    int numpoints;    /**< Number of points defining the body.            */
    double **coord;  /** points' coordinates. */
    double  s [3];    /**< Support mapping computed last.                 */
};

/**
 * @brief Structure for a simplex.
 */
struct simplex {
    int    nvrtx   ;       /**< Number of simplex's vertices. 			*/
    double vrtx    [4][3]; /**< Coordinates of simplex's vertices. 		*/
    int    wids    [4];    /**< Label of the simplex's vertices. 			*/
    double lambdas [4];    /**< Barycentric coordiantes for each vertex.  */
    double  p [4][3];     /**< Points of P that form the simplex  */
    double  q [4][3];	    /**< Points of Q that form the simplex  */
};


/**
 * @brief The GJK algorithm which returns the minimum distance between
 * two bodies.
 */
extern double gjk( struct bd, struct bd, struct simplex * ) ;
