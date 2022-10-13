//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "random_generator.hpp"
#include <chrono>


RandomGenerator::RandomGenerator(int seed)
{
    vslNewStream(&_stream, VSL_BRNG_MT19937, seed);

    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, _stream, 1000, _numbers.data(), 0.0f, 1.0f);

    _index = 0;
}

RandomGenerator::~RandomGenerator()
{
    vslDeleteStream(&_stream);
}

float
RandomGenerator::get_float(float a,
                           float b)
{
    const auto u = _numbers[_index];

    _index++;

    _update();

    return a + (b - a) * u;
}

std::vector<float>
RandomGenerator::get_floats(float a,
                            float b,
                            int32_t c)
{
    std::vector<float> randoms{};

    for (int32_t i = 0; i < c; i++)
    {
        randoms.push_back(get_float(a, b));
    }
    return randoms;
}

int
RandomGenerator::get_integer(const int a,
                             const int b)
{
    int r[0];

    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, _stream, 1, r, a, b);

    return r[0];
}

std::vector<int>
RandomGenerator::get_integers(const int a,
                              const int b,
                              const int length)
{
    std::vector<int> r(length, 0);

    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, _stream, length, r.data(), a, b);

    return r;
}

void
RandomGenerator::_update()
{
    if (_index == 1000)
    {
        vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, _stream, 1000, _numbers.data(), 0.0f, 1.0f);

        _index = 0;
    }
}


std::ostream &
operator<<(std::ostream &output, const RandomGenerator &source)
{
    output << std::endl;

    output << "[RandomGenerator (Instance): " << &source << "]" << std::endl;

    return output;
}



