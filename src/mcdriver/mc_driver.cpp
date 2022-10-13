//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include <random>
#include <iostream>
#include <chrono>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <omp.h>

#include "mkl.h"
#include "random_generator.hpp"
#include "math_operations.hpp"
#include "mc_driver.hpp"
#include "struct_generator.hpp"


MCDriver::MCDriver()
{
}

MCDriver::~MCDriver()
{
}

void
MCDriver::compute(int32_t NBOX, const float BOX_SIZE, const float DENSITY_FACTOR , const int ITERATIONS,
                  const float DOSE_VALUE, const int SPLIT_SIZE, const float DOSE_THRESHOLD,
                  const float DOSE_REDUCTION,const float CONNECTING_DISTANCE, const float TRANSLATING_DISTANCE)
{

    for (int32_t b = 0; b < NBOX; b++)
    {
        RandomGenerator gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        auto t_start = std::chrono::high_resolution_clock::now();

        StructureGenerator sga(86.09f, 5.0f);
        auto box = sga.generate(std::to_string(b), DENSITY_FACTOR, 0.0f, 0.0f, 0.0f,
                                BOX_SIZE, BOX_SIZE, BOX_SIZE, 2.0, 5.0, 5.0);

        //std::vector<float> dose_values = read_dose(dose_file);
        box.irradiate(DOSE_VALUE, SPLIT_SIZE);

        float max_dose = box.get_irradiation_max_value();

        /*  statistics  */
        std::vector<double> timing;
        int initial_monomers, initial_polymers;
        int activations = 0, connections = 0, divisions = 0;
        int contry = 0, divtry = 0;
        initial_monomers = box.number_of_monomers();
        initial_polymers = box.number_of_polymers();
        std::vector<std::vector<int>> constats(51, std::vector<int>(3));
        std::vector<int> divstats(3);
        std::vector<int> rotstats(2);

        /*  cycle   */
        for (int32_t i = 0; i < ITERATIONS; i++)
        {
            float capture_coeff = _calculate_capture();
            int32_t rpoli, rmoni;
            float dose_value_at;
            if (check_activation_v1B(box, dose_value_at, rpoli, rmoni, DOSE_THRESHOLD, DOSE_REDUCTION))
            {
                activations++;

                float P = _calculate_damage(box.polymer(rpoli).number_of_monomers(), max_dose, dose_value_at);

                if (capture_coeff >= P)
                {
                    contry++;
                    if (apply_connection_rotation_translation(box, rpoli, rmoni, TRANSLATING_DISTANCE,
                                                              CONNECTING_DISTANCE, constats, rotstats))
                    {
                        connections++;
                    }
                } else
                {
                    divtry++;
                    if (apply_division_rotation_translation(box, rpoli, rmoni, TRANSLATING_DISTANCE, divstats,
                                                            rotstats))
                    {
                        divisions++;
                    }
                }
            }

            if (i % 5000 == 0 && i > 1)
            {
                auto t_end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> t_elapsed_gpu_final = t_end - t_start;

            prepare_statistics(box, i, initial_monomers, initial_polymers,
                               activations, contry, connections, divtry, divisions,
                               timing, constats, divstats, rotstats, t_elapsed_gpu_final);
        }
    }
}

}


bool
MCDriver::check_activation_v1B(UnitBox &box, float &dose_value, int32_t &polymer_index, int32_t &monomer_index,
                               float DOSE_THRESHOLD, float DOSE_REDUCTION)
{
    RandomGenerator gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto dose_loc_index = gen.get_integer(0, box.get_irradiation().size() - 1);
    float dosex, dosey, dosez;
    box.get_irradiation(dose_loc_index, dosex, dosey, dosez, dose_value);

    if (box.find_nearest_for_radical(dosex, dosey, dosez, 12.0f, polymer_index, monomer_index))
    {
        if (dose_value > DOSE_THRESHOLD)
        {
            box.change_irradiation(dose_loc_index, dose_value - DOSE_REDUCTION);
            return true;
        }
    }
    return false;
}


float
MCDriver::_calculate_damage(int32_t mon_length, float max_dose, float cell_dose)
{
    float U = cell_dose / max_dose;
    float P;

    switch (mon_length)
    {
        case 1:
            P = 1.00f * ((1.0f - 0.0f) + 0.0f * std::erf(0.00f * U));
            break;
        case 2:
            P = 0.83f * ((1.0f - 0.8f) + 0.8f * std::erf(1.00f * U));
            break;
        case 3:
            P = 0.72f * ((1.0f - 0.8f) + 0.8f * std::erf(5.20f * U));
            break;
        case 4:
            P = 0.59f * ((1.0f - 0.4f) + 0.4f * std::erf(3.71f * U));
            break;
        case 5:
            P = 0.51f * ((1.0f - 0.6f) + 0.6f * std::erf(1.28f * U));
            break;
        case 6 ... 10:
            P = 0.23f * ((1.0f - 0.75f) + 0.75f * std::erf(0.42f * U * U));
            break;
        case 11 ... 15:
            P = 0.2f * ((1.0f - 0.67f) + 0.67f * std::erf(0.36f * U * U));
            break;
        default:
            P = 0.18f * ((1.0f - 0.6f) + 0.6f * std::erf(0.33f * U * U));
    }
    return P;
}

float
MCDriver::_calculate_capture()
{
    RandomGenerator gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    return (gen.get_float(0.0f, 1.0f));
}

bool
MCDriver::apply_rotation_translation(UnitBox &box, int32_t rpoli, float transl_dist)
{
    Polymer activated_polymer = box.polymer(rpoli);
    Branched_polymer activated_branpol, original_branpol;

    if (activated_polymer.is_branched())
    {
        activated_branpol = box.find_branpol(activated_polymer.get_bid());
        original_branpol = activated_branpol;
        box.remove_branpol(original_branpol);
        Rotransl rotransl;
        rotransl.rotate_translate(activated_branpol, transl_dist);
        if (box.check_branpol_intersection(activated_branpol))
        {
            box.add_branpol(original_branpol);
            return false;
        } else
        {
            box.add_branpol(activated_branpol);
            return true;
        }
    } else   /* for simple polymer */
    {
        Polymer original_polymer(activated_polymer);
        box.remove_polymer(original_polymer);
        Rotransl rotransl;
        rotransl.rotate_translate(activated_polymer, transl_dist);
        if (box.check_polymer_intersection(activated_polymer) || box.check_boundaries(activated_polymer))
        {
            box.add(original_polymer);
            return false;
        } else
        {
            box.add(activated_polymer);
            return true;
        }
    }
}

bool
MCDriver::internal_rotation_translation_branpol(UnitBox &box, Branched_polymer &activated_branpol, float transl_dist)
{
    Rotransl rotransl;
    rotransl.rotate_translate(activated_branpol, transl_dist);
    if (box.check_branpol_intersection(activated_branpol))
    {
        return false;
    } else
    {
        return true;
    }
}


bool
MCDriver::internal_rotation_translation_polymer(UnitBox &box, Polymer &activated_polymer, float transl_dist)
{
    Rotransl rotransl;
    rotransl.rotate_translate(activated_polymer, transl_dist);
    if (box.check_polymer_intersection(activated_polymer) || box.check_boundaries(activated_polymer))
    {
        return false;
    } else
    {
        return true;
    }
}


bool
MCDriver::apply_division_rotation_translation(UnitBox &box, int32_t rpoli, int32_t rmoni, float transl_dist,
                                              std::vector<int> &divstats, std::vector<int> &rotstats)
{
    Polymer activated_polymer = box.polymer(rpoli);
    Branched_polymer activated_branpol, original_branpol;
    divstats[0]++;

    if (activated_polymer.number_of_monomers() == 1 || rmoni == 0)
    {
        rotstats[0]++;
        if (apply_rotation_translation(box, rpoli, transl_dist))
        {
            rotstats[1]++;
            return false;
        }
        return false;
    }

    if (activated_polymer.is_branched())
    {

        activated_branpol = box.find_branpol(activated_polymer.get_bid());
        original_branpol = activated_branpol;
        box.remove_branpol(original_branpol);

        Division division;
        std::vector<Branched_polymer> b_output;
        std::vector<Polymer> p_output;
        if (division.perform_branpol_division_rotransl(activated_branpol, activated_polymer, rmoni, transl_dist,
                                                       b_output, p_output))
        {

            auto check_group = _get_polymers_for_check(b_output, p_output);
            if (box.check_group_intersection(check_group))
            {

                box.add_branpol(original_branpol);
                return false;
            } else
            {
                box.remove_branpol(original_branpol);
                for (auto &b: b_output) box.add_branpol(b);
                for (auto &p: p_output) box.add(p);
                divstats[2]++;
                return true;
            }
        } else
        {

            Branched_polymer rot_branpol(original_branpol);
            rotstats[0]++;
            if (internal_rotation_translation_branpol(box, rot_branpol, transl_dist))
            {
                rotstats[1]++;
                box.add_branpol(rot_branpol);
                return false;
            } else
            {
                box.add_branpol(original_branpol);
                return false;
            }
        }
    } else   /* for simple polymer */
    {

        Polymer original_polymer(activated_polymer);
        box.remove_polymer(original_polymer);

        Division division;
        std::vector<Polymer> output;
        if (division.perform_polymer_division_rotransl(activated_polymer, rmoni, transl_dist, output))
        {

            if (!box.check_polymer_intersection(output[0]) && !box.check_polymer_intersection(output[1]) &&
                !box.check_boundaries(output[0]) && !box.check_boundaries(output[1]))
            {

                box.add(output[0]);
                box.add(output[1]);
                divstats[1]++;
                return true;
            } else
            {

                box.add(original_polymer);
                return false;
            }
        } else
        {
            Polymer rotpol(original_polymer);
            rotstats[0]++;
            if (internal_rotation_translation_polymer(box, rotpol, transl_dist))
            {
                rotstats[1]++;
                box.add(rotpol);
                return false;
            }
            box.add(original_polymer);
            return false;
        }
    }
}


bool
MCDriver::apply_connection_rotation_translation(UnitBox &box, int32_t rpoli, int32_t rmoni, float transl_dist,
                                                float connect_distance,
                                                std::vector<std::vector<int>> &constats, std::vector<int> &rotstats)
{

    int32_t tpoli, tmoni;
    if (box.find_nearest_skip_same_branch(rpoli, rmoni, connect_distance, tpoli, tmoni))
    {
        Polymer act_pol = box.polymer(rpoli);
        Polymer targ_pol = box.polymer(tpoli);

        std::string option;
        Polymer act_pol_cpy, targ_pol_cpy;
        Branched_polymer act_branpol, act_branpol_cpy, targ_branpol, targ_branpol_cpy;
        if (act_pol.is_branched())
        {
            option = "B";
            act_branpol = box.find_branpol(act_pol.get_bid());
            act_branpol_cpy = act_branpol;
            box.remove_branpol(act_branpol);
        } else
        {
            option = "N";
            act_pol_cpy = act_pol;
            box.remove_polymer(act_pol);
        }
        if (targ_pol.is_branched())
        {
            option += "B";
            targ_branpol = box.find_branpol(targ_pol.get_bid());
            targ_branpol_cpy = targ_branpol;
            box.remove_branpol(targ_branpol);
        } else
        {
            option += "N";
            targ_pol_cpy = targ_pol;
            box.remove_polymer(targ_pol);
        }
        Connection conn;
        bool intersection = false;
        std::vector<Polymer> pol_output;
        std::vector<Branched_polymer> branpol_output;
        if (conn.perform_scenario(box, constats, act_pol, rmoni, act_branpol, targ_pol, tmoni, targ_branpol,
                                  pol_output, branpol_output))
        {
            Rotransl rot;
            for (auto &p: pol_output)
            {
                if (!intersection)
                {
                    rot.rotate_translate(p, transl_dist);
                    if (box.check_polymer_intersection(p) || box.check_boundaries(p)) intersection = true;
                }
            }
            for (auto &b: branpol_output)
            {
                if (!intersection)
                {
                    rot.rotate_translate(b, transl_dist);
                    if (box.check_branpol_intersection(b)) intersection = true;
                }
            }
            if (!intersection)
            {
                for (auto &p: pol_output) box.add(p);
                for (auto &b: branpol_output) box.add_branpol(b);
                return true;
            }
        }
        if (option == "NN")
        {
            Polymer rot_pol = act_pol_cpy;
            rotstats[0]++;
            if (internal_rotation_translation_polymer(box, rot_pol, transl_dist))
            {
                rotstats[1]++;
                box.add(rot_pol);
                box.add(targ_pol_cpy);
            } else
            {
                box.add(act_pol_cpy);
                box.add(targ_pol_cpy);
            }
        } else if (option == "BN")
        {
            Branched_polymer rot_branpol = act_branpol_cpy;
            rotstats[0]++;
            if (internal_rotation_translation_branpol(box, rot_branpol, transl_dist))
            {
                rotstats[1]++;
                box.add_branpol(rot_branpol);
                box.add(targ_pol_cpy);
            } else
            {
                box.add_branpol(act_branpol_cpy);
                box.add(targ_pol_cpy);
            }
        } else if (option == "NB")
        {
            Polymer rot_pol = act_pol_cpy;
            rotstats[0]++;
            if (internal_rotation_translation_polymer(box, rot_pol, transl_dist))
            {
                rotstats[1]++;
                box.add(rot_pol);
                box.add_branpol(targ_branpol_cpy);
            } else
            {
                box.add(act_pol_cpy);
                box.add_branpol(targ_branpol_cpy);
            }
        } else if (option == "BB")
        {
            Branched_polymer rot_branpol = act_branpol_cpy;
            rotstats[0]++;
            if (internal_rotation_translation_branpol(box, rot_branpol, transl_dist))
            {
                rotstats[1]++;
                box.add_branpol(rot_branpol);
                box.add_branpol(targ_branpol_cpy);
            } else
            {
                box.add_branpol(act_branpol_cpy);
                box.add_branpol(targ_branpol_cpy);
            }
        }
        return false;
    } else
    {
        rotstats[0]++;
        if (apply_rotation_translation(box, rpoli, transl_dist))
        {
            rotstats[1]++;
            return false;
        }
    }
    return false;
}


void
MCDriver::prepare_statistics(UnitBox &box, int box_nr, int initial_monomers, int initial_polymers, int activations,
                             int contry, int connections, int divtry, int divisions,
                             std::vector<double> &timing, std::vector<std::vector<int>> &constats,
                             std::vector<int> &divstats, std::vector<int> &rostats,
                             std::chrono::duration<double> t_elapsed_gpu_final)
{
    std::cout << "Preparing statistics" << std::endl;
    box.save("result_box_" + std::to_string(box_nr) + ".txt");
    std::ofstream my_file("statistics_" + std::to_string(box_nr) + ".txt");
    my_file << t_elapsed_gpu_final.count()
            << std::endl;
    my_file << "Initial monomers: " << initial_monomers << " Result monomers: " << box.number_of_monomers()
            << std::endl;
    my_file << "Initial_polymers : " << initial_polymers << " Result polymers: " << box.number_of_polymers()
            << std::endl;
    float sum_lengths = 0;
    for (auto p = 0; p < box.number_of_polymers(); p++)
        sum_lengths += box.polymer(p).number_of_monomers();
    my_file << "Average length of result polymers: " << sum_lengths / box.number_of_polymers() << std::endl;
    my_file << "Activations: " << activations << std::endl;
    my_file << "Connections (incl rot): attempts " << contry << " Succeeded " << connections << std::endl;
    my_file << "Divisions (incl rot): attempts " << divtry << " Succeeded " << divisions << std::endl;
    my_file << "Rotations: attempts " << rostats[0] << " Succeeded " << rostats[1] << std::endl;
    std::cout << "Preparing branching statistics" << std::endl;
    my_file << "Branching statistics calculation " << std::endl;
    std::vector<float> average_number, average_length, correlation;
    box.calculate_branching_alternative(average_number, average_length, correlation);
    my_file << "Average branch number: " << std::endl;
    for (auto &i: average_number) my_file << i << "  ";
    my_file << std::endl;
    my_file << "Average branch length: " << std::endl;
    for (auto &i: average_length) my_file << i << "  ";
    my_file << std::endl;
    my_file << "Correlation: " << std::endl;
    for (auto &i: average_number) my_file << i << "  ";
    my_file << std::endl;
    std::cout << "Preparing connections statistics" << std::endl;
    my_file << "Detail Connections: " << std::endl;
    for (auto &s: constats) my_file << s[0] << ":" << s[1] << "-" << s[2] << ", ";
    my_file << std::endl;
    std::cout << "Preparing division statistics" << std::endl;
    my_file << "Detail Divisions: " << std::endl;
    for (auto &s: divstats) my_file << s << " : ";
    my_file << std::endl;
}


int32_t
MCDriver::_find_element_in_vec(Polymer &polymer, std::vector<Polymer> &polymers)
{
    auto it = std::find(polymers.begin(), polymers.end(), polymer);
    if (it != polymers.end())
        return (it - polymers.begin());
    else
        return (-1);
}

std::vector<Polymer>
MCDriver::_get_polymers_for_check(std::vector<Branched_polymer> &branpols, std::vector<Polymer> &pols)
{
    std::vector<Polymer> group;
    for (auto &p: pols) group.push_back(p);
    for (auto &b: branpols)
    {
        auto branches = b.get_branches();
        auto parent = b.get_parent();
        group.insert(group.end(), branches.begin(), branches.end());
        group.push_back(parent);
    }
    return group;
}
