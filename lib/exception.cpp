// Custom exception types.
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

#include "exception.h"

invalid_file_format::invalid_file_format(const std::string& what_arg)
    : what_arg_(what_arg.c_str())
{
}

invalid_file_format::invalid_file_format(const char* what_arg)
    : what_arg_(what_arg)
{
}

const char* invalid_file_format::what() const noexcept
{
    return what_arg_;
}

unsupported_file_format::unsupported_file_format(const std::string& what_arg)
    : what_arg_(what_arg.c_str())
{
}

unsupported_file_format::unsupported_file_format(const char* what_arg)
    : what_arg_(what_arg)
{
}

const char* unsupported_file_format::what() const noexcept
{
    return what_arg_;
}
