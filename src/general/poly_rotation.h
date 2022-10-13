//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

/* helper class to rotate the polymer structures */

using Mat = std::vector<std::vector<double> >;

struct Pt
{
    double x, y, z;

    Pt(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z)
    {}
};

inline Pt operator+(const Pt &a, const Pt &b)
{ return {a.x + b.x, a.y + b.y, a.z + b.z}; }

inline Pt operator-(const Pt &a, const Pt &b)
{ return {a.x - b.x, a.y - b.y, a.z - b.z}; }

inline Pt operator/(const Pt &a, double d)
{ return {a.x / d, a.y / d, a.z / d}; }

inline Pt operator*(double d, const Pt &a)
{ return {d * a.x, d * a.y, d * a.z}; }

inline Pt cross(const Pt &a, const Pt &b)
{ return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }

inline double dot(const Pt &a, const Pt &b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

inline double normsq(const Pt &a)
{ return dot(a, a); }

inline double len(const Pt &a)
{ return sqrt(normsq(a)); }

inline std::ostream &operator<<(std::ostream &out, const Pt a)
{ return out << a.x << " " << a.y << " " << a.z << " "; }

inline Pt operator*(const Mat &M, const Pt &v)
{
    return {M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z,
            M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z,
            M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z};
}

inline
Mat operator*(const Mat &M, const Mat &N)
{
    Mat R(M.size(), std::vector<double>(N[0].size(), 0.0));
    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < N[0].size(); j++)
        {
            for (int k = 0; k < N.size(); k++) R[i][j] += M[i][k] * N[k][j];
        }
    }
    return R;
}


inline
Mat transpose(const Mat &M)
{
    Mat R(M.size(), std::vector<double>(M[0].size(), 0.0));
    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < M[0].size(); j++)
        {
            R[j][i] = M[i][j];
        }
    }
    return R;
}


inline
Mat rotationMatrix(const Pt &axis, float radians)
{
    Mat result(3, std::vector<double>(3, 0.0));
    Pt n = axis / len(axis);
    double cosA = cos(radians), sinA = sin(radians), cosA1 = 1.0 - cosA;
    result[0][0] = cosA + n.x * n.x * cosA1;
    result[0][1] = +n.x * n.y * cosA1 - n.z * sinA;
    result[0][2] = +n.x * n.z * cosA1 + n.y * sinA;
    result[1][0] = +n.y * n.x * cosA1 + n.z * sinA;
    result[1][1] = cosA + n.y * n.y * cosA1;
    result[1][2] = +n.y * n.z * cosA1 - n.x * sinA;
    result[2][0] = +n.z * n.x * cosA1 - n.y * sinA;
    result[2][1] = +n.z * n.y * cosA1 + n.x * sinA;
    result[2][2] = cosA + n.z * n.z * cosA1;
    return result;
}

inline
void axisAngle(const Mat &R, Pt &axis, float &angle)
{
    float angle_p = (0.5f * (R[0][0] + R[1][1] + R[2][2] - 1));
    if (angle_p > 1) angle_p = 1.00f;
    if (angle_p < -1) angle_p = -1.00f;
    angle = acos(angle_p);
    float s = sin(angle);

    if (fabs(s) < 1.0e-6)
    {
        int i = 0;
        if (R[1][1] > R[i][i]) i = 1;
        if (R[2][2] > R[i][i]) i = 2;
        axis = Pt(R[0][i], R[1][i], R[2][i]);
        if (i == 0) axis.x += 1;
        if (i == 1) axis.y += 1;
        if (i == 2) axis.z += 1;
    } else
    {
        axis = (0.5 / s) * Pt(R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]);
    }
    axis = axis / len(axis);
}