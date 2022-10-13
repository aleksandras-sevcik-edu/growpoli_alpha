//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include <chrono>
#include "connection.hpp"

bool Connection::perform_scenario(UnitBox &box, std::vector<std::vector<int>> &con_stats,
                                  Polymer &activated, int32_t act_moni, Branched_polymer &act_branpol,
                                  Polymer &target, int32_t targ_moni, Branched_polymer &targ_branpol,
                                  std::vector<Polymer> &pol_output, std::vector<Branched_polymer> &branpol_output)
{
    for (auto i = 0; i < 51; i++) con_stats[i][0] = i;
    con_stats[50][1]++;

    /* calculating which scenario */
    std::string scenario;
    auto anom = activated.number_of_monomers();
    auto tnom = target.number_of_monomers();
    scenario += (activated.is_branched()) ? "B" : "M";
    scenario += (anom < 10) ? "S" : "L";
    if (anom >= 10)
    {
        scenario += ((act_moni < 4) || (anom - act_moni < 4)) ? "E" : "M";
    }
    scenario += "-";
    scenario += (target.is_branched()) ? "B" : "M";
    scenario += (tnom < 10) ? "S" : "L";
    if (tnom >= 10)
    {
        if (target.is_branched())
        {
            auto parent = targ_branpol.get_parent();
            if (target == parent)
                scenario += "P";
            else
                scenario += "M";
        } else if ((targ_moni < 4) || ((tnom - targ_moni) < 4))
            scenario += "E";
        else
            scenario += "M";
    }

    if (scenario == "MS-MS")                                    /* scenario No 1 */
    {
        con_stats[1][1]++;
        if (merge(box, activated, target, false))
        {
            pol_output.push_back(activated);

            con_stats[1][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MLE-MS")                                    /* scenario No 2 */
    {
        con_stats[2][1]++;
        if (merge(box, activated, target, false))
        {
            pol_output.push_back(activated);

            con_stats[2][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MLM-MS")                                    /* scenario No 3 */
    {
        con_stats[3][1]++;
        if (mono_branch(box, activated, act_moni, target, targ_moni))
        {
            Branched_polymer bp(activated, target);
            branpol_output.push_back(bp);

            con_stats[3][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BS-MS")                                    /* scenario No 4 */
    {
        con_stats[4][1]++;
        auto reserve(activated);
        if (merge(box, activated, target, true))
        {
            if (act_branpol.change(reserve, activated))
            {
                branpol_output.push_back(act_branpol);

                con_stats[4][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLE-MS")                                    /* scenario No 5 */
    {
        con_stats[5][1]++;
        auto reserve(activated);
        if (merge(box, activated, target, true))
        {
            if (act_branpol.change(reserve, activated))
            {
                branpol_output.push_back(act_branpol);

                con_stats[5][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLM-MS")                                    /* scenario No 6 */
    {
        con_stats[6][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        if (mono_branch(box, activated, act_moni, target, targ_moni))
        {
            if (act_branpol.add_branch(activated, act_moni, target))
            {
                branpol_output.push_back(act_branpol);

                con_stats[6][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLP-MS")                                    /* scenario No 7 */
    {
        con_stats[7][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        if (mono_branch(box, activated, act_moni, target, targ_moni))
        {
            if (act_branpol.add_branch(activated, act_moni, target))
            {
                branpol_output.push_back(act_branpol);
                con_stats[7][2]++;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MS-MLE")                                    /* scenario No 8 */
    {
        con_stats[8][1]++;
        if (merge(box, target, activated, false))
        {
            pol_output.push_back(target);
            con_stats[8][2]++;
            return true;
        }
        return false;
    }
    if (scenario == "MLE-MLE")                                    /* scenario No 9 */
    {
        con_stats[9][1]++;
        if (merge(box, activated, target, false))
        {
            pol_output.push_back(activated);
            con_stats[9][2]++;
            return true;
        }
        return false;
    }
    if (scenario == "MLM-MLE")                                    /* scenario No 10 */
    {
        con_stats[10][1]++;
        if (mono_branch(box, activated, act_moni, target, targ_moni))
        {
            Branched_polymer bp(activated, target);
            branpol_output.push_back(bp);
            con_stats[10][2]++;
            return true;
        }
        return false;
    }
    if (scenario == "BS-MLE")                                    /* scenario No 11 */
    {
        con_stats[11][1]++;
        auto reserve(activated);
        if (merge(box, activated, target, true))
        {
            if (act_branpol.change(reserve, activated))
            {
                branpol_output.push_back(act_branpol);
                con_stats[11][2]++;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLE-MLE")                                    /* scenario No 12 */
    {
        con_stats[12][1]++;
        auto reserve(activated);
        if (merge(box, activated, target, true))
        {
            if (act_branpol.change(reserve, activated))
            {
                branpol_output.push_back(act_branpol);
                con_stats[12][2]++;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLM-MLE")                                    /* scenario No 13 */
    {
        con_stats[13][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        if (mono_branch(box, activated, act_moni, target, targ_moni))
        {
            if (act_branpol.add_branch(activated, act_moni, target))
            {
                branpol_output.push_back(act_branpol);
                con_stats[13][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLP-MLE")                                    /* scenario No 14 */
    {
        con_stats[14][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        if (mono_branch(box, activated, act_moni, target, targ_moni))
        {
            if (act_branpol.add_branch(activated, act_moni, target))
            {
                branpol_output.push_back(act_branpol);

                con_stats[14][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MS-MLM")                                    /* scenario No 15 */
    {
        con_stats[15][1]++;
        if (mono_branch(box, target, targ_moni, activated, targ_moni))
        {
            Branched_polymer bp(target, activated);
            branpol_output.push_back(bp);
            con_stats[15][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MLE-MLM")                                    /* scenario No 16 */
    {
        con_stats[16][1]++;
        if (mono_branch(box, target, targ_moni, activated, targ_moni))
        {
            Branched_polymer bp(target, activated);
            branpol_output.push_back(bp);
            con_stats[16][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MLM-MLM")                                    /* scenario No 17 */
    {
        con_stats[17][1]++;
        std::vector<Polymer> result;
        if (double_branch(box, activated, act_moni, target, targ_moni, result))
        {
            Branched_polymer bp(activated, result);
            branpol_output.push_back(bp);
            con_stats[17][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BS-MLM")                                    /* scenario No 18 */
    {
        con_stats[18][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, activated, act_moni, target, targ_moni, result))
        {
            if (act_branpol.add_branches(activated, act_moni, result))
            {
                branpol_output.push_back(act_branpol);
                con_stats[18][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLE-MLM")                                    /* scenario No 19 */
    {
        con_stats[19][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, activated, act_moni, target, targ_moni, result))
        {
            if (act_branpol.add_branches(activated, act_moni, result))
            {
                branpol_output.push_back(act_branpol);
                con_stats[19][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLM-MLM")                                    /* scenario No 20 */
    {
        con_stats[20][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, activated, act_moni, target, targ_moni, result))
        {
            if (act_branpol.add_branches(activated, act_moni, result))
            {
                branpol_output.push_back(act_branpol);
                con_stats[20][2]++;
                std::cout << scenario << " Ok!" << std::endl;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BLP-MLM")                                    /* scenario No 21 */
    {
        con_stats[21][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, activated, act_moni, target, targ_moni, result))
        {
            if (act_branpol.add_branches(activated, act_moni, result))
            {
                branpol_output.push_back(act_branpol);
                con_stats[21][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MS-BS")                                    /* scenario No 22 */
    {
        con_stats[22][1]++;
        auto reserve(target);
        if (merge(box, target, activated, true))
        {
            if (targ_branpol.change(reserve, target))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[22][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLE-BS")                                    /* scenario No 23 */
    {
        con_stats[23][1]++;
        auto reserve(target);
        if (merge(box, target, activated, true))
        {
            if (targ_branpol.change(reserve, target))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[23][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLM-BS")                                    /* scenario No 24 */
    {
        con_stats[24][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, target, targ_moni, activated, act_moni, result))
        {
            if (targ_branpol.add_branches(target, targ_moni, result))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[24][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BS-BS")                                    /* scenario No 25 */
    {
        con_stats[25][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[25][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLE-BS")                                    /* scenario No 26 */
    {
        con_stats[26][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[26][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLM-BS")                                    /* scenario No 27 */
    {
        con_stats[27][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Branched_polymer> bresults;
        std::vector<Polymer> presults;
        if (branch_swap(box, activated, act_moni, act_branpol, target, targ_moni, targ_branpol,
                        bresults, presults))
        {
            branpol_output.push_back(bresults[0]);
            branpol_output.push_back(bresults[1]);
            con_stats[27][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLP-BS")                                    /* scenario No 28 */
    {
        con_stats[28][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[28][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MS-BLE")                                    /* scenario No 29 */
    {
        con_stats[29][1]++;
        auto reserve(target);
        if (merge(box, target, activated, true))
        {
            if (targ_branpol.change(reserve, target))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[29][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLE-BLT")                                    /* scenario No 30 */
    {
        con_stats[30][1]++;
        auto reserve(target);
        if (merge(box, target, activated, true))
        {
            if (targ_branpol.change(reserve, target))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[30][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLM-BLE")                                    /* scenario No 31 */
    {
        con_stats[31][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, target, targ_moni, activated, act_moni, result))
        {
            if (targ_branpol.add_branches(target, targ_moni, result))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[31][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BS-BLE")                                    /* scenario No 32 */
    {
        con_stats[32][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, target, targ_moni, targ_branpol,
                         activated, act_moni, act_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[32][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLE-BLE")                                    /* scenario No 33 */
    {
        con_stats[33][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[33][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLM-BLE")                                    /* scenario No 34 */
    {
        con_stats[34][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Branched_polymer> bresults;
        std::vector<Polymer> presults;
        if (branch_swap(box, activated, act_moni, act_branpol, target, targ_moni, targ_branpol,
                        bresults, presults))
        {
            branpol_output.push_back(bresults[0]);
            branpol_output.push_back(bresults[1]);
            con_stats[34][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLP-BLE")                                    /* scenario No 35 */
    {
        con_stats[35][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[35][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MS-BLM")                                    /* scenario No 36 */
    {
        con_stats[36][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        if (mono_branch(box, target, targ_moni, activated, targ_moni))
        {
            if (targ_branpol.add_branch(target, targ_moni, activated))
            {
                branpol_output.push_back(targ_branpol);

                con_stats[36][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLE-BLM")                                    /* scenario No 37 */
    {
        con_stats[37][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        if (mono_branch(box, target, targ_moni, activated, targ_moni))
        {
            if (targ_branpol.add_branch(target, targ_moni, activated))
            {
                branpol_output.push_back(targ_branpol);

                con_stats[37][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLM-BLM")                                    /* scenario No 38 */
    {
        con_stats[38][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, target, targ_moni, activated, act_moni, result))
        {
            if (targ_branpol.add_branches(target, targ_moni, result))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[38][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BS-BLM")                                    /* scenario No 39 */
    {
        con_stats[39][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Branched_polymer> bresults;
        std::vector<Polymer> presults;
        if (branch_swap(box, activated, act_moni, act_branpol, target, targ_moni, targ_branpol,
                        bresults, presults))
        {
            branpol_output.push_back(bresults[0]);
            branpol_output.push_back(bresults[1]);
            con_stats[39][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLE-BLM")                                    /* scenario No 40 */
    {
        con_stats[40][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Branched_polymer> bresults;
        std::vector<Polymer> presults;
        if (branch_swap(box, activated, act_moni, act_branpol, target, targ_moni, targ_branpol,
                        bresults, presults))
        {
            branpol_output.push_back(bresults[0]);
            branpol_output.push_back(bresults[1]);
            con_stats[40][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLM-BLM")                                    /* scenario No 41 */
    {
        con_stats[41][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        std::vector<Branched_polymer> bresults;
        std::vector<Polymer> presults;
        if (branch_swap(box, activated, act_moni, act_branpol, target, targ_moni, targ_branpol,
                        bresults, presults))
        {
            branpol_output.push_back(bresults[0]);
            branpol_output.push_back(bresults[1]);
            con_stats[41][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLP-BLM")                                    /* scenario No 42 */
    {
        con_stats[42][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[42][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "MS-BLP")                                    /* scenario No 43 */
    {
        con_stats[43][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        if (mono_branch(box, target, targ_moni, activated, targ_moni))
        {
            if (targ_branpol.add_branch(target, targ_moni, activated))
            {
                branpol_output.push_back(targ_branpol);

                con_stats[43][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLE-BLP")                                    /* scenario No 44 */
    {
        con_stats[44][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        if (mono_branch(box, target, targ_moni, activated, targ_moni))
        {
            if (targ_branpol.add_branch(target, targ_moni, activated))
            {
                branpol_output.push_back(targ_branpol);

                con_stats[44][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "MLM-BLP")                                    /* scenario No 45 */
    {
        con_stats[45][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        std::vector<Polymer> result;
        if (double_branch(box, target, targ_moni, activated, act_moni, result))
        {
            if (targ_branpol.add_branches(target, targ_moni, result))
            {
                branpol_output.push_back(targ_branpol);
                con_stats[45][2]++;
                con_stats[50][2]++;;
                return true;
            }
        }
        return false;
    }
    if (scenario == "BS-BLP")                                    /* scenario No 46 */
    {
        con_stats[46][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, target, targ_moni, targ_branpol,
                         activated, act_moni, act_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[46][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLE-BLP")                                    /* scenario No 47 */
    {
        con_stats[47][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, target, targ_moni, targ_branpol,
                         activated, act_moni, act_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[47][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLM-BLP")                                    /* scenario No 48 */
    {
        con_stats[48][1]++;
        if (targ_branpol.is_monomer_occupied(target, targ_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, target, targ_moni, targ_branpol,
                         activated, act_moni, act_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[48][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    if (scenario == "BLP-BLP")                                    /* scenario No 49 */
    {
        con_stats[49][1]++;
        if (act_branpol.is_monomer_occupied(activated, act_moni)) return false;
        Branched_polymer result;
        if (branch_merge(box, activated, act_moni, act_branpol,
                         target, targ_moni, targ_branpol, result))
        {
            branpol_output.push_back(result);
            con_stats[49][2]++;
            con_stats[50][2]++;;
            return true;
        }
        return false;
    }
    abort();
}


bool Connection::merge(UnitBox &box, Polymer &activated, Polymer &target, bool is_branched)
{
    Polymer basePolymer(activated);
    Polymer nearestPolymer(target);

    float px, py, pz, pxl, pyl, pzl, ox, oy, oz, oxl, oyl, ozl;
    activated.position_center(px, py, pz, 0);
    activated.position_center(pxl, pyl, pzl, activated.number_of_monomers() - 1);
    target.position_center(ox, oy, oz, 0);
    target.position_center(oxl, oyl, ozl, target.number_of_monomers() - 1);
    float begin_begin = (((ox - px) * (ox - px)) + ((oy - py) * (oy - py)) + ((oz - pz) * (oz - pz)));
    float begin_end = (((oxl - px) * (oxl - px)) + ((oyl - py) * (oyl - py)) + ((ozl - pz) * (ozl - pz)));
    float end_end = (((oxl - pxl) * (oxl - pxl)) + ((oyl - pyl) * (oyl - pyl)) + ((ozl - pzl) * (ozl - pzl)));
    float end_begin = (((ox - pxl) * (ox - pxl)) + ((oy - pyl) * (oy - pyl)) + ((oz - pzl) * (oz - pzl)));
    std::vector<float> pos = {begin_begin, begin_end, end_begin, end_end};
    auto it = std::min_element(std::begin(pos), std::end(pos));
    int index = std::distance(std::begin(pos), it);

    if (!is_branched)
    {
        if (index == 0)
        {
            basePolymer.reverse();
        } else if (index == 1)
        {
            basePolymer.reverse();
            nearestPolymer.reverse();
        } else if (index == 2)
        {
            nearestPolymer.reverse();
        }
    }

    ExponentialMap endchain = basePolymer.rotation(basePolymer.number_of_monomers() - 1);
    ExponentialMap begchain = nearestPolymer.rotation(0);
    ExponentialMap join = _align_bound(endchain, begchain, 1.22173f);
    float fx, fy, fz;
    nearestPolymer.position(fx, fy, fz, 0);
    nearestPolymer.rotate_polymer({fx, fy, fz}, join);
    float rx, ry, rz;
    basePolymer.position(rx, ry, rz, basePolymer.number_of_monomers());
    nearestPolymer.set_origins(rx - nearestPolymer.monomer(0).length() / 2,
                               ry - nearestPolymer.monomer(0).width() / 2, rz);
    basePolymer.polymer_merge(nearestPolymer);
    basePolymer.calculate_coords();

    if (basePolymer.intersecting_itself())
    {
        return false;
    }

    activated = basePolymer;
    return true;
}


bool Connection::mono_branch(UnitBox &box, Polymer &activated, int32_t mon_index, Polymer &target, int32_t targ_moni)
{
    if (mon_index == 0 || mon_index == (activated.number_of_monomers() - 1)) return false;

    if (targ_moni >= target.number_of_monomers() - 4) target.reverse();


    ExponentialMap endchain = activated.rotation(mon_index);
    ExponentialMap begchain = target.rotation(0);
    ExponentialMap join = _align_bound(endchain, begchain, 4.36332);

    float fx, fy, fz;
    target.position(fx, fy, fz, 0);
    target.rotate_polymer({fx, fy, fz}, join);
    float rx, ry, rz;
    activated.position(rx, ry, rz, mon_index);
    target.set_origins(rx - target.monomer(0).length() / 2,
                       ry - target.monomer(0).width() / 2, rz);

    target.calculate_coords();
    if (!_is_intersecting_from_second_mon(activated, target)) return true;

    return false;
}


bool Connection::double_branch(UnitBox &box, Polymer &activated, int32_t act_mon, Polymer &target, int32_t targ_mon,
                               std::vector<Polymer> &result)
{
    if (act_mon == 0 || act_mon == (activated.number_of_monomers() - 1)) return false;
    if (targ_mon == 0 || targ_mon == (target.number_of_monomers() - 1)) return false;
    if (target.number_of_monomers() < 2) return false;

    Polymer first_half, second_half;
    target.divide(targ_mon, first_half, second_half);
    first_half.reverse();

    ExponentialMap endchain = activated.rotation(act_mon);
    ExponentialMap begchain_f = first_half.rotation(first_half.number_of_monomers() - 1);
    ExponentialMap begchain_s = second_half.rotation(0);
    ExponentialMap join_s = _align_bound(endchain, begchain_s, 4.36332);
    ExponentialMap join_f = _align_bound(endchain, begchain_f, 4.36332);

    float sx, sy, sz;
    second_half.position(sx, sy, sz, 0);
    second_half.rotate_polymer({sx, sy, sz}, join_s);
    float fx, fy, fz;
    first_half.position(fx, fy, fz, 0);
    first_half.rotate_polymer({fx, fy, fz}, join_f);

    float rx, ry, rz;
    activated.position(rx, ry, rz, act_mon);
    second_half.set_origins(rx - second_half.monomer(0).length() / 2,
                            ry - second_half.monomer(0).width() / 2, rz);

    float rx2, ry2, rz2;
    activated.position(rx2, ry2, rz2, act_mon + 1);
    first_half.set_origins(rx2 - first_half.monomer(0).length() / 2,
                           ry2 - first_half.monomer(0).width() / 2, rz2);

    first_half.calculate_coords();
    second_half.calculate_coords();
    if (_is_intersecting_from_second_mon(activated, first_half)) return false;
    if (_is_intersecting_from_second_mon(activated, second_half)) return false;
    if (_is_intersecting_from_first_mon(first_half, second_half)) return false;

    result.push_back(first_half);
    result.push_back(second_half);

    return true;
}


bool Connection::branch_merge(UnitBox &box, Polymer &activated, int32_t act_mon, Branched_polymer &act_branpol,
                              Polymer &target, int32_t targ_mon, Branched_polymer &targ_branpol,
                              Branched_polymer &result)
{
    if ((act_mon == 0) || (act_mon == activated.number_of_monomers() - 1)) return false;

    int junction = target.number_of_monomers() - 1;

    std::cout << target.number_of_monomers() << " " << targ_mon << std::endl;
    ExponentialMap endchain = activated.rotation(act_mon);
    ExponentialMap begchain = target.rotation(junction);
    ExponentialMap join = _align_bound(endchain, begchain, 4.36332);
    float rx, ry, rz;
    target.position(rx, ry, rz, junction);
    auto target_label = target.get_label();
    auto tpolymers = targ_branpol.to_polymers();
    for (auto &pol: tpolymers) pol.rotate_polymer({rx, ry, rz}, join);
    auto rot_target = _find_polymer(target_label, tpolymers);
    //move
    float ax, ay, az;
    activated.position(ax, ay, az, act_mon);
    _move_polymers(rot_target, junction + 1, ax, ay, az, tpolymers);

    auto apolymers = act_branpol.to_polymers();
    apolymers.insert(apolymers.end(), tpolymers.begin(), tpolymers.end());
    for (auto &p: apolymers) p.calculate_coords();
    result.from_polymers(apolymers);
    result.make_longest_parent();

    if (result.is_intersecting_between()) return false;

    return true;
}


bool Connection::branch_swap(UnitBox &box, Polymer &activated, int32_t act_mon, Branched_polymer &act_branpol,
                             Polymer &target, int32_t targ_mon, Branched_polymer &targ_branpol,
                             std::vector<Branched_polymer> &bresults, std::vector<Polymer> &presults)
{
    //rotations
    ExponentialMap endchain_1 = activated.rotation(act_mon);
    ExponentialMap begchain_1 = target.rotation(targ_mon);
    ExponentialMap join_1 = _align_bound(endchain_1, begchain_1, 4.36332);
    ExponentialMap endchain_2 = target.rotation(targ_mon);
    ExponentialMap begchain_2 = activated.rotation(act_mon);
    ExponentialMap join_2 = _align_bound(endchain_2, begchain_2, 4.36332);
    //positions
    float ax, ay, az, tx, ty, tz;
    activated.position(ax, ay, az, act_mon);
    target.position(tx, ty, tz, targ_mon);
    //division
    std::vector<Branched_polymer> first_div_bran, sec_div_bran;
    std::vector<Polymer> first_div_pol, sec_div_pol, first_specials, sec_specials;
    std::string scen;
    bool fd = act_branpol.make_division_for_swap(activated, act_mon, first_div_bran, first_div_pol, first_specials);
    bool sd = targ_branpol.make_division_for_swap(target, targ_mon, sec_div_bran, sec_div_pol, sec_specials);

    if (first_div_bran.size() == 2)
    {
        scen += "BB";
    } else if (first_div_bran.size() == 1 && first_div_pol.size() == 1)
    {
        scen += "BP";
    }
    if (sec_div_bran.size() == 2)
    {
        scen += "BB";
    } else if (sec_div_bran.size() == 1 && sec_div_pol.size() == 1)
    {
        scen += "BP";
    }
    Branched_polymer result_a, result_b;
    std::vector<Polymer> first_1, first_2, second_1, second_2;
    std::string label_1, label_2;

    if (fd && sd && scen == "BBBB")
    {
        first_1 = first_div_bran[0].to_polymers();
        first_2 = first_div_bran[1].to_polymers();
        second_1 = sec_div_bran[0].to_polymers();
        second_2 = sec_div_bran[1].to_polymers();
        label_1 = first_specials[1].get_label();
        label_2 = sec_specials[1].get_label();
    } else if (fd && sd && scen == "BBBP")
    {
        first_1 = first_div_bran[0].to_polymers();
        first_2 = first_div_bran[1].to_polymers();
        second_1 = sec_div_bran[0].to_polymers();
        second_2 = sec_div_pol;
        label_1 = first_specials[1].get_label();
        label_2 = sec_div_pol[0].get_label();
    } else if (fd && sd && scen == "BPBB")
    {
        first_1 = first_div_bran[0].to_polymers();
        first_2 = first_div_pol;
        second_1 = sec_div_bran[0].to_polymers();
        second_2 = sec_div_bran[1].to_polymers();
        label_1 = first_div_pol[0].get_label();
        label_2 = sec_specials[1].get_label();
    } else if (fd && sd && scen == "BPBP")
    {
        first_1 = first_div_bran[0].to_polymers();
        first_2 = first_div_pol;
        second_1 = sec_div_bran[0].to_polymers();
        second_2 = sec_div_pol;
        label_1 = first_div_pol[0].get_label();
        label_2 = sec_div_pol[0].get_label();
    } else
        return false;

    //rotate
    for (auto &pol: first_2) pol.rotate_polymer({ax, ay, az}, join_2);
    for (auto &pol: second_2) pol.rotate_polymer({tx, ty, tz}, join_1);
    //move
    auto first_2_act = _find_polymer(label_1, first_2);
    _move_polymers(first_2_act, 0, tx, ty, tz, first_2);
    auto second_2_act = _find_polymer(label_2, second_2);
    _move_polymers(second_2_act, 0, ax, ay, az, second_2);


    first_1.insert(first_1.end(), second_2.begin(), second_2.end());
    second_1.insert(second_1.end(), first_2.begin(), first_2.end());

    for (auto &p: first_1) p.calculate_coords();
    result_a.from_polymers(first_1);
    for (auto &p: second_1) p.calculate_coords();
    result_b.from_polymers(second_1);

    result_a.make_longest_parent();
    result_b.make_longest_parent();

    if (result_a.is_intersecting_between()) return false;
    if (result_b.is_intersecting_between()) return false;

    bresults.push_back(result_a);
    bresults.push_back(result_b);

    return true;
}


bool
Connection::_is_intersecting_from_second_mon(Polymer &activated, Polymer &target)
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();
    if (activated.polymer_intersecting_from_second(target, &bd1, &bd2))
    {
        free_bd(&bd1);
        free_bd(&bd2);
        return true;
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}


bool
Connection::_is_intersecting_from_first_mon(Polymer &activated, Polymer &target)
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();
    if (activated.polymer_intersecting(target, &bd1, &bd2))
    {
        free_bd(&bd1);
        free_bd(&bd2);
        return true;
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}

void
Connection::_move_polymers(Polymer &target, int32_t act_mon, float fx, float fy, float fz,
                           std::vector<Polymer> &polymers)
{
    //initial position
    float ix, iy, iz, ixo, iyo, izo, nixo, niyo, nizo;
    target.position(ixo, iyo, izo, 0);
    target.position(ix, iy, iz, act_mon);
    target.set_origins(fx - target.monomer(0).length() / 2,
                       fy - target.monomer(0).width() / 2, fz);
    target.position(nixo, niyo, nizo, 0);
    float difx = nixo - ixo - (ix - ixo);
    float dify = niyo - iyo - (iy - iyo);
    float difz = nizo - izo - (iz - izo);
    for (auto &polymer: polymers)
    {
        float cx, cy, cz;
        polymer.position(cx, cy, cz, 0);
        polymer.set_origins(cx - polymer.monomer(0).length() / 2 + difx,
                            cy - polymer.monomer(0).width() / 2 + dify,
                            cz + difz);
    }
}


ExponentialMap
Connection::_align_bound(ExponentialMap &previous_rot, ExponentialMap &align_to, float delta) const
{
    auto prng = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    //const float delta = 1.22173f; //( 180 - 110 ) * PI / 180.0 = 1.22173
    const float pi = 3.14159265358979323846f;
    Pt ex = {1, 0, 0}, ey = {0, 1, 0}, ez = {0, 0, 1};
    float px, py, pz, ptheta;
    previous_rot.get_rotation(px, py, pz, ptheta);
    Pt previous = {px, py, pz};
    Mat R = rotationMatrix(previous, ptheta);
    previous = R * ez;
    Pt perp = cross(previous, ex);
    if (len(perp) < 0.01f) perp = cross(previous, ey);
    R = rotationMatrix(previous, prng.get_float(0.0f, 2 * pi));
    perp = R * perp;
    previous = previous / len(previous);
    perp = perp / len(perp);
    Pt axis = cos(delta) * previous + sin(delta) * perp;
    Pt edge = sin(delta) * previous - cos(delta) * perp;
    R = rotationMatrix(axis, prng.get_float(0.0f, 2 * pi));
    edge = R * edge;
    Pt third = cross(axis, edge);
    Mat A = {{edge.x, third.x, axis.x},
             {edge.y, third.y, axis.y},
             {edge.z, third.z, axis.z}};
    // alignment
    float ax, ay, az, atheta;
    align_to.get_rotation(ax, ay, az, atheta);
    Mat B = rotationMatrix({ax, ay, az}, atheta);
    Mat Binv = transpose(B);
    Mat ABinv = A * Binv;
    Pt monomer_axis;
    float new_theta;
    axisAngle(ABinv, monomer_axis, new_theta);

    return (ExponentialMap(monomer_axis.x, monomer_axis.y, monomer_axis.z, new_theta));
}

bool
Connection::_is_centers_too_close(Polymer &activated, int32_t amoni, Polymer &target)
{
    float ax, ay, az, tx, ty, tz;
    activated.position_center(ax, ay, az, amoni);
    target.position_center(tx, ty, tz, 0);

    if (((ax - tx) * (ax - tx) + (ay - ty) * (ay - ty) + (az - tz) * (az - tz)) <
        activated.monomer(amoni).width() * activated.monomer(amoni).width())
        return true;

    return false;
}


int
Connection::_find_position(const Polymer &polymer, std::vector<Polymer> &polymers)
{
    auto iter = std::find(polymers.begin(), polymers.end(), polymer);
    return (std::distance(polymers.begin(), iter));
}


Polymer
Connection::_find_polymer(const std::string &label, std::vector<Polymer> &polymers)
{
    for (auto &p: polymers)
    {
        if (p.get_label() == label)
            return p;
    }
    Polymer empty;
    return empty;
}
