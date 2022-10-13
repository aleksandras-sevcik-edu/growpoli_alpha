//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "struct_generator.hpp"
#include <iostream>
#include <cmath>
#include <chrono>
#include "math_operations.hpp"
#include <random>
#include <omp.h>


StructureGenerator::StructureGenerator(const float &monomer_mass,
                                       const float &monomer_concentration)

        : _monomer_mass(monomer_mass), _monomer_concentration(monomer_concentration)
{
}

StructureGenerator::~StructureGenerator()
{
}

UnitBox
StructureGenerator::generate(const std::string &label_box,
                             float density_factor,
                             const float box_origin_x,
                             const float box_origin_y,
                             const float box_origin_z,
                             const float box_x,
                             const float box_y,
                             const float box_z,
                             const float mlen,
                             const float mwid,
                             const float mhig)
{
    std::vector<Polymer> polymers;

    auto punits = _get_number_of_monomer_units(_monomer_mass, _monomer_concentration, density_factor,
                                               box_x - box_origin_x, box_y - box_origin_y,
                                               box_z - box_origin_z);
    while (punits > 0)
    {
        auto clen = _get_polymer_chain_length();

        clen = (clen >= punits) ? punits : clen;

        auto npol = _get_full_polymer(clen, mlen, mwid, mhig,
                                      box_origin_x, box_x,
                                      box_origin_y, box_y,
                                      box_origin_z, box_z);
        if (polymers.empty())
        {
            polymers.push_back(npol);
            punits -= clen;

        } else
        {
            if (!_check_polymer_intersection(polymers, npol))
            {
                polymers.push_back(npol);
                punits -= clen;
            }
        }
    }
    return UnitBox(polymers, {}, label_box, box_origin_x, box_origin_y, box_origin_z, box_x, box_y, box_z, {});
}


Polymer
StructureGenerator::_get_full_polymer(const float clen,
                                      const float mlen,
                                      const float mwid,
                                      const float mhig,
                                      const float box_origin_x,
                                      const float box_x,
                                      const float box_origin_y,
                                      const float box_y,
                                      const float box_origin_z,
                                      const float box_z) const
{
    while (true)
    {
        auto counter = 0;

        auto polymer = _connect_polymer(clen, mlen, mwid, mhig);

        while (true)
        {
            counter += 1;

            if (counter > 200)
            {
                break;
            }

            auto prng = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());

            float x_origin = prng.get_float(box_origin_x, box_x);

            float y_origin = prng.get_float(box_origin_y, box_y);

            float z_origin = prng.get_float(box_origin_z, box_z);

            polymer.set_origins(x_origin, y_origin, z_origin);

            bool is_out = false;

            for (int32_t i = 0; i < polymer.number_of_monomers(); i++)
            {
                float iposx, iposy, iposz;

                polymer.position(iposx, iposy, iposz, i);

                if ((iposx >= box_x) || (iposy >= box_y) || (iposz >= box_z)
                    || (iposx <= box_origin_x) || (iposy <= box_origin_y) || (iposz <= box_origin_z))
                {
                    is_out = true;

                    for (int32_t l = 0; l < polymer.number_of_monomers(); l++) polymer.fill_coords(l);

                    break;
                }
            }
            if (is_out)
            {
                continue;

            } else
            {
                for (int32_t l = 0; l < polymer.number_of_monomers(); l++) polymer.fill_coords(l);

                return polymer;
            }
        }
    }
}


Polymer
StructureGenerator::_connect_polymer(int32_t clen,
                                     float mlen,
                                     float mwid,
                                     float mhig) const
{
    auto out_checker = 0;

    auto checker = 1;

    int32_t size = 50;

    int32_t count = clen / size;

    int32_t size_last = clen - size * count;

    size = (size >= clen) ? clen : size;

    auto polymer_base = _get_polymer(mlen, mwid, mhig, size);

    if (polymer_base.number_of_monomers() == clen) return polymer_base;

    while (checker < count)
    {
        float posx, posy, posz;

        polymer_base.position(posx, posy, posz, polymer_base.number_of_monomers());

        auto add_polymer = _get_polymer(mlen, mwid, mhig, size);

        add_polymer.set_origins(posx - mlen / 2, posy - mwid / 2, posz);

        for (int32_t i = 0; i < add_polymer.number_of_monomers(); i++) add_polymer.fill_coords(i);

        if (!polymer_base.polymer_intersecting_connect(add_polymer))
        {
            polymer_base.polymer_merge(add_polymer);

            checker += 1;

            out_checker = 0;
        } else
        {
            out_checker += 1;
        }
        if (out_checker > 20)
        {
            if (polymer_base.number_of_monomers() > (2 * size))
            {
                polymer_base.shorten(2 * size);

                checker = checker - 2;

                if (checker == 0) checker = 1;

                out_checker = 0;
            } else
            {
                polymer_base = _get_polymer(mlen, mwid, mhig, size);

                checker = 1;

                out_checker = 0;
            }
        }
    }
    if (size_last > 0)
    {
        while (true)
        {
            float posx, posy, posz;

            polymer_base.position(posx, posy, posz, polymer_base.number_of_monomers());

            auto add_polymer = _get_polymer(mlen, mwid, mhig, size_last);

            add_polymer.set_origins(posx - mlen / 2, posy - mwid / 2, posz);

            for (int32_t i = 0; i < add_polymer.number_of_monomers(); i++) add_polymer.fill_coords(i);

            if (!polymer_base.polymer_intersecting_connect(add_polymer))
            {
                polymer_base.polymer_merge(add_polymer);

                break;
            }
        }
    }
    return polymer_base;
}


bool
StructureGenerator::_check_polymer_intersection(std::vector<Polymer> &polymers,
                                                const Polymer &polymer)
{
    auto ptr_polymer = &polymer;
    auto npol = polymers.size();
    auto ptr_polymers = polymers.data();

    const auto nthreads = omp_get_max_threads();

    bool check = false;

#pragma omp parallel shared (ptr_polymer, ptr_polymers, npol, check)
    {
#pragma omp single nowait
        {
            const auto batch_size = npol / nthreads;
            for (int32_t i = 0; i < nthreads; i++)
            {
                const auto bstart = batch_size * i;
                const auto bend = ((bstart + batch_size) > npol) ? npol : bstart + batch_size;

#pragma omp task firstprivate(i, bstart, bend) shared (check)
                {
                    struct bd bd1{}, bd2{};
                    bd1 = allocate_bd();
                    bd2 = allocate_bd();

                    for (auto j = bstart; j < bend; j++)
                    {
                        bool loc_check;
#pragma omp atomic read
                        loc_check = check;
                        if (loc_check) break;

                        if (ptr_polymer->polymer_intersecting(ptr_polymers[j], &bd1, &bd2))
                        {
#pragma omp atomic write
                            check = true;
                        }
                    }
                    free_bd(&bd1);
                    free_bd(&bd2);
                }
            }
        }
    }
    return check;
}


Polymer
StructureGenerator::_get_polymer(float mlen,
                                 float mwid,
                                 float mhig,
                                 int32_t nmon) const
{

    Polymer cpol;

    while (true)
    {
        Polymer npol;

        auto counter = 0;

        while (npol.number_of_monomers() < nmon)
        {
            auto numon = npol.number_of_monomers();

            Monomer mon({"m"}, Rectangle3D(mlen, mwid, mhig), {}, 0.0f, 0.0f, 0.0f);

            if (npol.number_of_monomers() == 0)
            {
                npol.add(mon, _random_rotation());
            } else
            {
                ExponentialMap prot = npol.rotation(npol.number_of_monomers() - 1);
                npol.add(mon, _random_rotation_bound(prot));
            }

            if (numon == npol.number_of_monomers()) counter += 1;

            if (numon + 1 == npol.number_of_monomers()) counter = 0;

            if (counter > 100) break;
        }

        if (npol.number_of_monomers() == nmon)
        {
            cpol = npol;

            cpol.fill_coords(cpol.number_of_monomers() - 1);

            break;
        }
    }
    return cpol;
}

ExponentialMap
StructureGenerator::_random_rotation() const
{
    const auto m_pi = getPiValue();
    auto prng = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    float x = prng.get_float(0.0f, 1.0f);
    float z = prng.get_float(0.0f, 1.0f);
    float y = prng.get_float(0.0f, 1.0f);

    auto r = ExponentialMap(0.0f, 0.0f, 1.0f, prng.get_float(0.0f, 2 * m_pi));

    float refx, refy, refz;

    r.rotate(refx, refy, refz, x, y, z);

    return ExponentialMap(refx, refy, refz, prng.get_float(0.0f, 2 * m_pi));
}


ExponentialMap
StructureGenerator::_random_rotation_bound(ExponentialMap &previous_rot) const
{
    auto prng = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    const float delta = 1.22173f; //( 180 - 110 ) * PI / 180.0 = 1.22173
    const float pi = getPiValue();

    Pt ex = {1, 0, 0}, ey = {0, 1, 0}, ez = {0, 0, 1};

    float px, py, pz, ptheta;
    previous_rot.get_rotation(px, py, pz, ptheta);
    Pt previous = {px, py, pz};                 // previous ROTATION axis
    // Convert it to vector which is the rotated z axis
    Mat R = rotationMatrix(previous, ptheta);
    previous = R * ez;                            // Should be the alignment of the monomer

    Pt perp = cross(previous, ex);
    if (len(perp) < 0.01) perp = cross(previous, ey);

    R = rotationMatrix(previous, prng.get_float(0.0f, 2 * pi));
    perp = R * perp;
    previous = previous / len(previous);  // normalization to length 1
    perp = perp / len(perp);              // ditto
    Pt axis = cos(delta) * previous + sin(delta) * perp;  // new rotated z axis

    Pt edge = sin(delta) * previous - cos(delta) * perp;
    R = rotationMatrix(axis, prng.get_float(0.0f, 2 * pi));
    edge = R * edge;

    Pt third = cross(axis, edge);

    /*        Rot(ex)  Rot(ey)  Rot(ez)
           (      |      |         |   )
       M = (    edge   third     axis  )
           (      |      |         |   )    */

    Mat M = {{edge.x, third.x, axis.x},
             {edge.y, third.y, axis.y},
             {edge.z, third.z, axis.z}};

    Pt monomer_axis;
    float new_theta;
    axisAngle(M, monomer_axis, new_theta);

    return ExponentialMap(monomer_axis.x, monomer_axis.y, monomer_axis.z, new_theta);
}


int32_t
StructureGenerator::_get_number_of_monomer_units(const float monomer_mass,
                                                 const float monomer_concentration,
                                                 const float density_factor,
                                                 const float box_x,
                                                 const float box_y,
                                                 const float box_z) const
{
    //return box_x * box_y * box_z * 1E-24 * monomer_concentration / 100 * density_factor * 6.02214076E23 / monomer_mass;
    return (box_x * box_y * box_z * 0.000001f * 47 * density_factor);
}


int32_t
StructureGenerator::_get_polymer_chain_length() const
{
    auto prng = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    float prob = prng.get_float(0.0f, 1.0f);

    if (prob < 0.33f)
        return 1;
    else if (prob > 0.63f)
        return 2;
    else
        return 3;
}


float
StructureGenerator::get_monomer_concentration() const
{
    return _monomer_concentration;
}


float
StructureGenerator::get_monomer_mass() const
{
    return _monomer_mass;
}


UnitBox
StructureGenerator::generate_test(const std::string &label_box,
                                  int punits,
                                  int clen,
                                  const float box_origin_x,
                                  const float box_origin_y,
                                  const float box_origin_z,
                                  const float box_x,
                                  const float box_y,
                                  const float box_z,
                                  const float mlen,
                                  const float mwid,
                                  const float mhig)
{
    std::vector<Polymer> polymers;

    while (punits > 0)
    {
        clen = (clen >= punits) ? punits : clen;

        auto npol = _get_full_polymer(clen, mlen, mwid, mhig,
                                      box_origin_x, box_x,
                                      box_origin_y, box_y,
                                      box_origin_z, box_z);
        if (polymers.empty())
        {
            polymers.push_back(npol);

            punits -= clen;

        } else
        {
            if (!_check_polymer_intersection(polymers, npol))
            {
                polymers.push_back(npol);

                punits -= clen;
            }
        }
    }

    return UnitBox(polymers, {}, label_box, box_origin_x, box_origin_y, box_origin_z, box_x, box_y, box_z, {});
}
