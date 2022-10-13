//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "rotransl.hpp"

bool Rotransl::rotate_translate(Polymer &act_polymer,
                                float translation_distance)
{
    std::vector<Polymer> polymers{act_polymer};

    if (act_polymer.number_of_monomers() > 1)
    {
        float mass_center_x, mass_center_y, mass_center_z;
        ExponentialMap first, second, third;
        _calculate_inertia(polymers, mass_center_x, mass_center_y, mass_center_z, first, second, third);

        for (auto &polymer: polymers)
        {
            polymer.rotate_polymer({mass_center_x, mass_center_y, mass_center_z}, first);
            polymer.rotate_polymer({mass_center_x, mass_center_y, mass_center_z}, second);
            polymer.rotate_polymer({mass_center_x, mass_center_y, mass_center_z}, third);
        }
    }
    RandomGenerator gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto randnr = gen.get_floats(-translation_distance, translation_distance, 3);
    for (auto &polymer: polymers)
    {
        float fx, fy, fz;
        polymer.position(fx, fy, fz, 0);
        polymer.set_origins(fx - polymer.monomer(0).length() / 2 + randnr[0],
                            fy - polymer.monomer(0).width() / 2 + randnr[1], fz + randnr[2]);
    }
    act_polymer = polymers[0];
    act_polymer.calculate_coords();

    return true;
}


bool Rotransl::rotate_translate(Branched_polymer &branpol,
                                float translation_distance)
{
    auto polymers = branpol.get_branches();
    polymers.push_back(branpol.get_parent());
    if (polymers.empty()) return false;

    float mass_center_x, mass_center_y, mass_center_z;
    ExponentialMap first, second, third;
    _calculate_inertia(polymers, mass_center_x, mass_center_y, mass_center_z, first, second, third);

    for (auto &polymer: polymers)
    {
        polymer.rotate_polymer({mass_center_x, mass_center_y, mass_center_z}, first);
        polymer.rotate_polymer({mass_center_x, mass_center_y, mass_center_z}, second);
        polymer.rotate_polymer({mass_center_x, mass_center_y, mass_center_z}, third);
    }

    RandomGenerator gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto randnr = gen.get_floats(-translation_distance, translation_distance, 3);
    for (auto &polymer: polymers)
    {
        float fx, fy, fz;
        polymer.position(fx, fy, fz, 0);
        polymer.set_origins(fx - polymer.monomer(0).length() / 2 + randnr[0],
                            fy - polymer.monomer(0).width() / 2 + randnr[1], fz + randnr[2]);
    }
    branpol.from_polymers(polymers);
    branpol.calculate_coords();

    return true;
}


void
Rotransl::_calculate_inertia(std::vector<Polymer> &polymers,
                             float &mass_center_x,
                             float &mass_center_y,
                             float &mass_center_z,
                             ExponentialMap &first_rot,
                             ExponentialMap &second_rot,
                             ExponentialMap &third_rot) const
{
    //compute mass center
    float sum_x = 0, sum_y = 0, sum_z = 0;
    int32_t nom = 0;
    for (auto &polymer: polymers)
    {
        for (int32_t i = 0; i < polymer.number_of_monomers(); i++)
        {
            float px, py, pz;
            polymer.position_center_alternative_2(px, py, pz, i);
            sum_x += px;
            sum_y += py;
            sum_z += pz;
        }
        nom += polymer.number_of_monomers();
    }
    auto nomd = 1 / nom;
    mass_center_x = sum_x * nomd;
    mass_center_y = sum_y * nomd;
    mass_center_z = sum_z * nomd;

    std::vector<float> axes, inertia;
    std::vector<float> tensor(9, 0.0);
    for (auto &polymer: polymers)
    {
        for (int32_t i = 0; i < polymer.number_of_monomers(); i++)
        {
            float xcoord, ycoord, zcoord;
            polymer.position_center_alternative_2(xcoord, ycoord, zcoord, i);
            const auto rx = xcoord - mass_center_x;
            const auto ry = ycoord - mass_center_y;
            const auto rz = zcoord - mass_center_z;
            const auto mt = 1.0f;
            tensor[0] += mt * (ry * ry + rz * rz);
            tensor[1] -= mt * rx * ry;
            tensor[2] -= mt * rx * rz;
            tensor[4] += mt * (rx * rx + rz * rz);
            tensor[5] -= mt * ry * rz;
            tensor[8] += mt * (rx * rx + ry * ry);
        }
    }
    tensor[3] = tensor[1];
    tensor[6] = tensor[2];
    tensor[7] = tensor[5];
    inertia = std::vector<float>(3, 0.0);
    axes = std::vector<float>(9, 0.0);
    std::vector<int32_t> idx(6);
    int32_t nval = 0, ndim = 3;
    LAPACKE_ssyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', ndim, tensor.data(),
                   ndim, 0.0, 0.0, 0, 0, 1.0e-13, &nval, inertia.data(),
                   axes.data(), ndim, idx.data());
    RandomGenerator gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::vector<ExponentialMap> rotations;
    for (int32_t j = 0; j < 3; j++)
    {
        if (inertia[j] <= 0) inertia[j] = 0.000001f;
        auto rf = gen.get_float(1.0f, std::log(inertia[j]) + inertia[j]);
        float teta = 360 * (1 - std::erf(inertia[j] / rf));
        ExponentialMap newr = {axes[j * 3], axes[j * 3 + 1], axes[j * 3 + 2], teta};
        rotations.push_back(newr);
    }
    first_rot = rotations[0];
    second_rot = rotations[1];
    third_rot = rotations[2];
}