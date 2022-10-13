//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include <cstdint>


#include "unit_box.hpp"
#include "random_generator.hpp"

/**
 Class StructureGenerator implements randomized unit box generator.
 */
class StructureGenerator
{
public:
    float _monomer_mass;

    float _monomer_concentration;

    int32_t _get_number_of_monomer_units(const float monomer_mass,
                                         const float monomer_concentration,
                                         const float density_factor,
                                         const float box_x,
                                         const float box_y,
                                         const float box_z) const;

    int32_t _get_polymer_chain_length() const;

    Polymer _get_polymer(float mlen,
                         float mwid,
                         float mhig,
                         int32_t nmon) const;

    ExponentialMap _random_rotation() const;

    ExponentialMap _random_rotation_bound(ExponentialMap &rot) const;

    Polymer _get_full_polymer(const float clen,
                              const float mlen,
                              const float mwid,
                              const float mhig,
                              const float box_origin_x,
                              const float box_x,
                              const float box_origin_y,
                              const float box_y,
                              const float box_origin_z,
                              const float box_z) const;

    Polymer _connect_polymer(int32_t clen,
                             float mlen,
                             float mwid,
                             float mhig) const;

    static bool _check_polymer_intersection(std::vector<Polymer> &polymers,
                                            const Polymer &polymer);

    StructureGenerator(const float &monomer_mass,
                       const float &monomer_concentration);

    ~StructureGenerator();

    UnitBox generate(const std::string &label_box,
                     const float density_factor,
                     const float box_origin_x,
                     const float box_origin_y,
                     const float box_origin_z,
                     const float box_x,
                     const float box_y,
                     const float box_z,
                     const float mlen,
                     const float mwid,
                     const float mhig);

    float get_monomer_concentration() const;

    float get_monomer_mass() const;

    bool check_junction(Polymer &polymer, Polymer &other) const;

    bool _check_polymer_intersection(const std::vector<Polymer> &polymers, const Polymer &polymer) const;

    UnitBox generate_test(const std::string &label_box, int punits, int clen, const float box_origin_x,
                          const float box_origin_y,
                          const float box_origin_z, const float box_x, const float box_y, const float box_z,
                          const float mlen,
                          const float mwid, const float mhig);
};
