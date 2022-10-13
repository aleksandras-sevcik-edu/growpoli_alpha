//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include "polymer.hpp"

/* Implements coordinate calculations for polymer structures */

namespace coord
{

    bool polymer_intersecting(const float *coords, const int *index, const int *size,
                              const float *pcoord, const int *pindex, const int *psize,
                              int jpol, struct bd *bd1, struct bd *bd2);

    void get_mcoord(std::vector<float> &coords, int index, struct bd *bd);

    void get_all_coords(const std::vector<Polymer> &polymers,
                        std::vector<float> &coords);

    int get_all_size(const std::vector<Polymer> &polymers, std::vector<int> &index, std::vector<int> &size);

}
