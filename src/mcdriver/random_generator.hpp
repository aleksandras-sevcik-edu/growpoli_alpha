//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once


#include <array>
#include <vector>
#include <ostream>
#include <cstdio>
#include "mkl_vsl.h"

class RandomGenerator
{
    VSLStreamStatePtr _stream;

    std::array<float, 1000> _numbers;

    int32_t _index;

    void _update();

public:

    RandomGenerator(int seed);

    ~RandomGenerator();

    float
    get_float(float a, float b);

    std::vector<float>
    get_floats(float a, float b, int32_t c);

    int
    get_integer(const int a,
                const int b);

    std::vector<int> get_integers(const int a,
                                  const int b,
                                  const int length);

    friend std::ostream &operator<<(std::ostream &output,
                                    const RandomGenerator &source);


};
