// Quick and dirty logging object.
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

#include "log.h"

int Log::global_level_ = 0;

void Log::SetLevel(int level)
{
    global_level_ = level;
}

Log& endl(Log& log)
{
    if (log.level_ >= log.global_level_)
        std::endl(std::cout);
    return log;
}

Log& Log::operator<<(Log& (*pf)(Log&))
{
    return pf(*this);
}
