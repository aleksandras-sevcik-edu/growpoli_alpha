#include <iostream>
#include "mc_driver.hpp"

const int32_t NBOX = 1;
const float BOX_SIZE = 500;
const float DENSITY_FACTOR = 1.0f;
const int ITERATIONS = 10000000;
const float DOSE_VALUE = 100.0f;
const int SPLIT_SIZE = 5;
const float DOSE_THRESHOLD = 1.0f;
const float DOSE_REDUCTION = 1.0f;
const float CONNECTING_DISTANCE = 15.0f;
const float TRANSLATING_DISTANCE = 30.0f;

int main()
{
    MCDriver mcd;
    mcd.compute(NBOX, BOX_SIZE, DENSITY_FACTOR, ITERATIONS, DOSE_VALUE,
                SPLIT_SIZE, DOSE_THRESHOLD, DOSE_REDUCTION, CONNECTING_DISTANCE,
                TRANSLATING_DISTANCE);
}