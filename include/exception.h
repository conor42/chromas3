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

#pragma once

#include <exception>
#include <string>

class invalid_file_format : public std::exception {
public:
    explicit invalid_file_format(const std::string& what_arg);
    explicit invalid_file_format(const char* what_arg);
    virtual const char* what() const noexcept;

private:
    const char* what_arg_;
};