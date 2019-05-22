//Copyright (c) 2018 Ultimaker B.V.


#include "Point3.h" //The headers we're implementing.

namespace arachne
{

Point3 Point3::operator +(const Point3& p) const
{
    return Point3(x + p.x, y + p.y, z + p.z);
}

Point3 Point3::operator -(const Point3& p) const
{
    return Point3(x - p.x, y - p.y, z - p.z);
}

Point3 Point3::operator *(const Point3& p) const
{
    return Point3(x * p.x, y * p.y, z * p.z);
}

Point3 Point3::operator /(const Point3& p) const
{
    return Point3(x / p.x, y / p.y, z / p.z);
}

Point3& Point3::operator +=(const Point3& p)
{
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
}

Point3& Point3::operator -=(const Point3& p)
{
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
}

Point3& Point3::operator *=(const Point3& p)
{
    x *= p.x;
    y *= p.y;
    z *= p.z;
    return *this;
}

Point3& Point3::operator /=(const Point3& p)
{
    x /= p.x;
    y /= p.y;
    z /= p.z;
    return *this;
}

bool Point3::operator ==(const Point3& p) const
{
    return x == p.x && y == p.y && z == p.z;
}

bool Point3::operator !=(const Point3& p) const
{
    return x != p.x || y != p.y || z != p.z;
}

}
