//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

/* helper class to use open_gjk */

#pragma once

#include <cstdio>
#include <vector>
#include "open_gjk.h"


inline bool
check_gjk_intersection(struct bd *bd1, struct bd *bd2)
{
    double dd;

    struct simplex s{};

    s.nvrtx = 0;

    dd = gjk(*bd1, *bd2, &s);

    //printf("Distance between bodies %f\n", dd);

    if (dd < 0.001)  return true;

    return false;
}


inline struct bd
allocate_bd()
{
    struct bd bd{};

    bd.numpoints = 8;

    bd.coord = (double **) malloc(8 * sizeof(double *));

    for (int32_t i = 0; i < 8; i++)
    {
        bd.coord[i] = (double *) malloc(3 * sizeof(double));
    }
    return bd;
}


inline void
free_bd(struct bd* bd)
{
    for (int32_t i = 0; i < 8; i++)
    {
        free(bd->coord[i]);
    }
    free(bd->coord);
}
