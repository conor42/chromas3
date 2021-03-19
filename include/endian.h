// Big-endian byte order conversions.
// 
// Copyright 2021 Conor N. McCarthy
//
// This file is part of Chromas 3.
//
// Chromas 3 is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chromas 3 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Chromas 3. If not, see < https://www.gnu.org/licenses/>.

#pragma once

#include <cstdlib>

inline unsigned short system_endian(unsigned short val)
{
#ifndef BIG_ENDIAN
    return _byteswap_ushort(val);
#else
    return val;
#endif
}

inline unsigned long system_endian(unsigned long val)
{
#ifndef BIG_ENDIAN
    return _byteswap_ulong(val);
#else
    return val;
#endif
}

inline unsigned int system_endian(unsigned int val)
{
#ifndef BIG_ENDIAN
    return _byteswap_ulong(val);
#else
    return val;
#endif
}

inline unsigned short system_endian(short val)
{
#ifndef BIG_ENDIAN
    return _byteswap_ushort(val);
#else
    return val;
#endif
}

inline unsigned long system_endian(long val)
{
#ifndef BIG_ENDIAN
    return _byteswap_ulong(val);
#else
    return val;
#endif
}

inline unsigned int system_endian(int val)
{
#ifndef BIG_ENDIAN
    return _byteswap_ulong(val);
#else
    return val;
#endif
}

template<typename T>
T read_bigendian(const char* src, size_t size)
{
#ifndef BIG_ENDIAN
    ptrdiff_t i = size - 1;
    T value = src[i];
    for (--i; i >= 0; --i)
        value = (value << 8) | src[i];
    return value;
#else
    T value = src[0];
    for (size_t i = 1; i < size; ++i)
        value = (value << 8) + src[i];
    return value;
#endif
}