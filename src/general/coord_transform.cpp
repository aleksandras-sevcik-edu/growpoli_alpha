//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "coord_transform.hpp"

namespace coord
{
    void get_mcoord(const float *coords, int index, struct bd *bd)
    {
        for (int32_t i = 0; i < 8; i++)
        {
            for (int32_t j = 0; j < 3; j++)
            {
                bd->coord[i][j] = coords[index + i*3 + j];
            }
        }
    }


    bool polymer_intersecting(const float *coords, const int *index, const int *size,
                              const float *pcoord, const int *pindex, const int *psize,
                              int jpol, struct bd *bd1, struct bd *bd2)
    {
        for (int32_t i = 0; i < psize[0]; i++)
        {
            get_mcoord(pcoord, 24 * i, bd1);

            for (int32_t j = 0; j < size[jpol]; j++)
            {
                get_mcoord(coords, index[jpol] + 24 * i, bd2);
                if (check_gjk_intersection(bd1, bd2))
                {
                    return true;
                }
            }
        }
        return false;
    }


    void get_all_coords(const std::vector<Polymer> &polymers,
                        std::vector<float> &coords)
    {
        int ids = 0;
        for (auto &p: polymers)
        {
            p.get_coords_to_vec(coords.data(), ids);
            ids+= 24 * p.number_of_monomers();
        }
    }

    int get_all_size(const std::vector<Polymer> &polymers, std::vector<int> &index, std::vector<int> &size)
    {
        index.clear();
        size.clear();

        index.reserve(polymers.size());
        size.reserve(polymers.size());

        int nsize = 0;
        for (auto &p : polymers )
        {
            const auto nmo = p.number_of_monomers();
            size.push_back(nmo);
            index.push_back(nsize);
            nsize += nmo * 24;
        }
        return nsize;
    }
}