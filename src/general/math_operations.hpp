//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

/*helper class for math constants and functions */

#include <iostream>
#include <vector>
#include <array>

inline float
getPiValue()
{
    return 3.14159265358979323846;
}

/* Boltzmann constant in eV⋅K−1 */
inline float
getKValue()
{
    return 0.00008617333262;
}

//vector cross product ( 2->1 * 2->3)
inline std::array<float, 3>
get_cross_product(float x1, float y1, float z1,
                  float x2, float y2, float z2,
                  float x3, float y3, float z3)
{
    std::array<float, 3> arr{};

    float a1 = x1 - x2;
    float b1 = y1 - y2;
    float c1 = z1 - z2;
    float a2 = x3 - x2;
    float b2 = y3 - y2;
    float c2 = z3 - z2;

    arr[0] = b1 * c2 - c1 * b2;
    arr[1] = c1 * a2 - a1 * c2;
    arr[2] = a1 * b2 - b1 * a2;

    return arr;
}


/* Checks if two given vectors are collinear and returns the equation of plane */
inline std::array<float, 4>
get_eop(float x1, float y1, float z1,
        float x2, float y2, float z2,
        float x3, float y3, float z3)
{
    float a1 = x1 - x2;
    float b1 = y1 - y2;
    float c1 = z1 - z2;
    float a2 = x3 - x2;
    float b2 = y3 - y2;
    float c2 = z3 - z2;

    std::array<float, 4> eop{};

    eop[0] = b1 * c2 - c1 * b2;
    eop[1] = c1 * a2 - a1 * c2;
    eop[2] = a1 * b2 - b1 * a2;
    eop[3] = (- eop[0] * x1 - eop[1] * y1 - eop[2] * z1);

    return eop;
}



