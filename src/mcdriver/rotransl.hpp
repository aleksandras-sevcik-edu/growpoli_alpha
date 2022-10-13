//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once
/* Implements rotation translation event*/

#include <chrono>

#include "mkl.h"
#include "unit_box.hpp"
#include "random_generator.hpp"

class Rotransl
{
public:

    bool rotate_translate(Polymer &polymer, float translation_distance);

    bool rotate_translate(Branched_polymer &branpol, float translation_distance);

    void
    _calculate_inertia(std::vector<Polymer> &polymers, float &mass_center_x, float &mass_center_y, float &mass_center_z,
                       ExponentialMap &first_rot, ExponentialMap &second_rot, ExponentialMap &third_rot) const;
};

