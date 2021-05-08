// ab1 file reader.
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

#include <cinttypes>
#include <ctime>
#include <memory>
#include <string>
#include <type_traits>
#include "endian.h"

class Ab1File
{
public:
    template <typename T>
    class Iterator
    {
    public:
        Iterator() : element_(nullptr), size_(0) {}

        Iterator(const char* data, size_t size) : element_(data), size_(size) {}

        Iterator(const Iterator& rval) : element_(rval.element_), size_(rval.size_) {}

        Iterator& operator+=(ptrdiff_t rval) {
            element_ += rval * size_; return *this;
        }

        Iterator& operator-=(ptrdiff_t rval) {
            element_ -= rval * size_; return *this;
        }

        Iterator& operator++() {
            element_ += size_; return *this;
        }

        Iterator& operator--() {
            element_ -= size_; return *this;
        }

        Iterator operator++(int) const {
            Iterator tmp(*this);
            element_ += size_;
            return tmp;
        }
        Iterator operator--(int) const {
            Iterator tmp(*this);
            element_ -= size_;
            return tmp;
        }

        const T* operator->() const {
            return reinterpret_cast<const T*>(element_);
        }

        // Endian conversion prevents the return of a reference here. Only valid for fundamental types and structs,
        // but nothing else is needed.
        const T operator[](ptrdiff_t rval) const {
            if (std::is_fundamental<T>::value)
                return read_bigendian<T>(element_ + rval * size_, size_);
            T value = 0;
            memcpy(&value, element_ + rval * size_, size_);
            return value;
        }

        const T operator*() const {
            return operator[](0);
        }

        ptrdiff_t operator-(const Iterator& rval) const {
            return (element_ - rval.element_) / size_;
        }

        Iterator operator+(ptrdiff_t rval) const {
            return Iterator(element_ + rval * size_, size_);
        }

        Iterator operator-(ptrdiff_t rval) const {
            return Iterator(element_ - rval * size_, size_);
        }

        friend Iterator operator+(ptrdiff_t lval, const Iterator& rval) {
            return Iterator(lval * size_ + rval.element_, size_);
        }

        friend Iterator operator-(ptrdiff_t lval, const Iterator& rval) {
            return Iterator(rval.element_ - lval * size_, size_);
        }

        bool operator==(const Iterator& rval) const {
            return element_ == rval.element_;
        }
        bool operator!=(const Iterator& rval) const {
            return element_ != rval.element_;
        }

    private:
        const char* element_;
        size_t size_;
    };

    enum Type : int16_t
    {
        BYTE = 1,
        CHAR = 2,
        USHORT= 3,
        SHORT = 4,
        LONG = 5,
        FLOAT = 7,
        DOUBLE = 8,
        DATE = 10,
        TIME = 11,
        BOOL = 13,
        PSTRING = 18,
        CSTRING = 19
    };

#pragma pack(push,1)

    // The spec declares signed types for some structure members for which negative values would not make sense,
    // so unsigned types are used here to simplify sanity checks.
    struct DirEntry
    {
        uint32_t name;
        int32_t number;
        Type data_type;
        uint16_t element_len;
        uint32_t elements;
        uint32_t bytes;
        uint32_t data;
        int32_t handle; // unused
    };

    struct Header
    {
        uint32_t magic; // 'ABIF'
        uint16_t version;
        DirEntry dir_pos;
        uint16_t reserved[47];

        Header(const char* buffer, size_t pos) {
            memcpy(&magic, buffer + pos, sizeof(Header));
        }
    };

    struct Date
    {
        int16_t year;
        uint8_t month;
        uint8_t day;
    };

    struct Time
    {
        uint8_t hour;
        uint8_t minute;
        uint8_t second;
        uint8_t hsecond;
    };

#pragma pack(pop)

    enum class SearchResult
    {
        SUCCESS,
        NOT_FOUND,
        INCOMPATIBLE_TYPE
    };

    Ab1File(const char* path);

    Ab1File(const Ab1File&) = delete;

    Ab1File(Ab1File&&) = delete;

    Ab1File& operator=(const Ab1File&) = delete;

    Ab1File& operator=(Ab1File&&) = delete;

    template <typename T>
    SearchResult SearchTag(const char* tag, int32_t number, Iterator<T>& begin, Iterator<T>& end)
    {
        DirEntry entry;
        SearchResult found = SearchTag(tag, number, entry);
        if (found != SearchResult::SUCCESS)
            return found;

        const char* data_ptr = file_buffer_.get() + entry.data;
        const char* data_end = data_ptr + size_t(entry.elements) * entry.element_len;
        switch (entry.data_type) {
        case BYTE:
        case CHAR:
        case USHORT:
        case SHORT:
        case LONG:
        case FLOAT:
        case DOUBLE:
        case BOOL:
            if(!std::is_fundamental<T>::value)
                return SearchResult::INCOMPATIBLE_TYPE;
            begin = Iterator<T>(data_ptr, entry.element_len);
            end = Iterator<T>(data_end, entry.element_len);
            break;

        case DATE:
        case TIME:
            if (std::is_fundamental<T>::value)
                return SearchResult::INCOMPATIBLE_TYPE;
            begin = Iterator<T>(data_ptr, entry.element_len);
            end = Iterator<T>(data_end, entry.element_len);
            break;

        default:
            return SearchResult::INCOMPATIBLE_TYPE;
        }
        return SearchResult::SUCCESS;
    }

    // Search for a 4-letter + number AB1 data section, perform endian conversion if necessary, and copy the contents into converted.
    // May return NOT_FOUND.
    SearchResult SearchTag(const char* tag, int32_t number, DirEntry& converted);

    // Search for a string (Pascal or C format) by tag and number, and copy it to result.
    // May return NOT_FOUND or INCOMPATIBLE_TYPE.
    SearchResult SearchTag(const char* tag, int32_t number, std::string& result);

    // Search for an integer value (any size) by tag and number, and copy it to result.
    // May return NOT_FOUND or INCOMPATIBLE_TYPE.
    SearchResult SearchTag(const char* tag, int32_t number, int32_t& value);

    // Search for date and time values by tags and number, and copy it to datetime.
    // May return NOT_FOUND or INCOMPATIBLE_TYPE.
    SearchResult DateTime(const char* date_tag, const char* time_tag, int32_t number, time_t& datetime);

private:

    void FindDirectory();

    DirEntry* SearchTag(const char* tag, int32_t number);

    size_t DataPos(DirEntry* entry);

    std::unique_ptr<char[]> file_buffer_;
    size_t file_start_;
    size_t file_size_;
    DirEntry* dir_begin_;
    DirEntry* dir_src_;
    DirEntry* dir_end_;
};