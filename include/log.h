// Quick and dirty logging object. I'll make a better one later.
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

#include <iostream>

class Log
{
public:
    Log(int level)
        : level_(level) {}

    template<typename T>
    Log& operator<< (T val) {
        if (level_ >= global_level_)
            std::cout.operator<<(val);
        return *this;
    }

    Log& operator<< (const char* s) {
        if (level_ >= global_level_)
            std::operator<<(std::cout, s);
        return *this;
    }

    Log& operator<< (char c) {
        if (level_ >= global_level_)
            std::operator<<(std::cout, c);
        return *this;
    }

    Log& operator<<(Log& (*pf)(Log&));

    friend Log& endl(Log& log);

    static void SetLevel(int level);

private:
    int level_;
    static int global_level_;
};

