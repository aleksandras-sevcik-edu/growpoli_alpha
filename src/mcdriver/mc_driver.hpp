//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include "unit_box.hpp"
#include "random_generator.hpp"
#include "rotransl.hpp"
#include "connection.hpp"
#include "division.hpp"

#include <cstdint>

/* Class MCDriver implements Monte Carlo driver for polymer irradiation process */

class MCDriver
{
public:

    MCDriver();
    ~MCDriver();

    void compute(int32_t NBOX, const float BOX_SIZE, const float DENSITY_FACTOR, const int ITERATIONS,
                 const float DOSE_VALUE, const int SPLIT_SIZE, const float DOSE_THRESHOLD, const float DOSE_REDUCTION,
                 const float CONNECTING_DISTANCE, const float TRANSLATING_DISTANCE);

    //std::vector<float> read_dose(const std::string &dose_file);

    float _calculate_damage(int32_t mon_length, float max_dose, float cell_dose);

    float _calculate_capture();

    static int32_t _find_element_in_vec(Polymer &polymer, std::vector<Polymer> &polymers);

    bool apply_rotation_translation(UnitBox &box, int32_t rpoli, float transl_dist);


    bool apply_division_rotation_translation(UnitBox &box, int32_t rpoli, int32_t rmoni, float transl_dist,
                                             std::vector<int> &divstats, std::vector<int> &rotstats);

    std::vector<Polymer> _get_polymers_for_check(std::vector<Branched_polymer> &branpols, std::vector<Polymer> &pols);

    bool apply_connection_rotation_translation(UnitBox &box, int32_t rpoli, int32_t rmoni, float transl_dist,
                                               float connect_distance, std::vector<std::vector<int>> &constats,
                                               std::vector<int> &rotstats);

    void prepare_statistics(UnitBox &box, int box_nr, int initial_monomers, int initial_polymers, int activations,
                            int contry, int connections, int divtry, int divisions,
                            std::vector<double> &timing, std::vector<std::vector<int>> &constats,
                            std::vector<int> &divstats, std::vector<int> &rostats,
                            std::chrono::duration<double> t_elapsed_gpu_final);

    bool internal_rotation_translation_branpol(UnitBox &box, Branched_polymer &activated_branpol, float transl_dist);

    bool internal_rotation_translation_polymer(UnitBox &box, Polymer &activated_polymer, float transl_dist);

    bool check_activation_v1B(UnitBox &box, float &dose_value, int32_t &polymer_index, int32_t &monomer_index,
                              float DOSE_THRESHOLD, float DOSE_REDUCTION);
};
