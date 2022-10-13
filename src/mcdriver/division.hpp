//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

/* Implements break-down event */

#include "unit_box.hpp"
#include "rotransl.hpp"


class Division
{
public:

    bool
    perform_polymer_division_rotransl(Polymer &polymer,
                                      int32_t moni,
                                      float transl_dist,
                                      std::vector<Polymer> &output);


    bool perform_branpol_division_rotransl(Branched_polymer &branpol,
                                           Polymer &polymer,
                                           int32_t moni,
                                           float transl_dist,
                                           std::vector<Branched_polymer> &output,
                                           std::vector<Polymer> &output_p);

    bool is_intersecting(std::vector<Branched_polymer> &output, std::vector<Polymer> &output_p);
};
