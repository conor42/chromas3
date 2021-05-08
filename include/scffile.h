// SCF file reader.
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
#include <string>
#include "endian.h"
#include "lookuptables.h"

class ScfFile
{
public:
    template <typename T>
    class Iterator
    {
    public:
        Iterator() : element_(nullptr), size_(0), pitch_(0) {}

        Iterator(const char* data, size_t size, size_t pitch) : element_(data), size_(size), pitch_(pitch) {}

        Iterator(const Iterator& rval) : element_(rval.element_), size_(rval.size_), pitch_(rval.pitch_) {}

        Iterator& operator+=(ptrdiff_t rval) {
            element_ += rval * pitch_; return *this;
        }

        Iterator& operator-=(ptrdiff_t rval) {
            element_ -= rval * pitch_; return *this;
        }

        Iterator& operator++() {
            element_ += pitch_; return *this;
        }

        Iterator& operator--() {
            element_ -= pitch_; return *this;
        }

        Iterator operator++(int) const {
            Iterator tmp(*this);
            element_ += pitch_;
            return tmp;
        }
        Iterator operator--(int) const {
            Iterator tmp(*this);
            element_ -= pitch_;
            return tmp;
        }

        const T* operator->() const {
            return reinterpret_cast<const T*>(element_);
        }

        // Endian conversion prevents the return of a reference here. Only valid for fundamental types and structs,
        // but nothing else is needed.
        const T operator[](ptrdiff_t rval) const {
            return read_bigendian<T>(element_ + rval * pitch_, size_);
        }

        const T operator*() const {
            return operator[](0);
        }

        ptrdiff_t operator-(const Iterator& rval) const {
            return (element_ - rval.element_) / pitch_;
        }

        Iterator operator+(ptrdiff_t rval) const {
            return Iterator(element_ + rval * pitch_, size_, pitch_);
        }

        Iterator operator-(ptrdiff_t rval) const {
            return Iterator(element_ - rval * pitch_, size_, pitch_);
        }

        friend Iterator operator+(ptrdiff_t lval, const Iterator& rval) {
            return Iterator(lval * pitch_ + rval.element_, size_, pitch_);
        }

        friend Iterator operator-(ptrdiff_t lval, const Iterator& rval) {
            return Iterator(rval.element_ - lval * pitch_, size_, pitch_);
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
        size_t pitch_;
    };

    template <typename T>
    class TraceIterator
    {
    public:
        TraceIterator()
            : element_(nullptr), size_(0), pitch_(0), delta_(0), sample_(0), transform_(false), init_(false)
        {}

        TraceIterator(const char* data, size_t size, size_t pitch, bool transform)
            : element_(data), size_(size), pitch_(pitch), delta_(0), sample_(0), transform_(transform), init_(true)
        {}

        TraceIterator(const TraceIterator& rval)
            : element_(rval.element_), size_(rval.size_), pitch_(rval.pitch_), delta_(rval.delta_), sample_(rval.sample_), transform_(rval.transform_), init_(rval.init_) {}

        TraceIterator& operator++() {
            if (init_)
                (void)operator*();
            element_ += pitch_;
            init_ = true;
            return *this;
        }

        TraceIterator operator++(int) {
            Iterator tmp(*this);
            operator++();
            return tmp;
        }

        const T operator*() {
            if (init_) {
                if (transform_) {
                    delta_ += read_bigendian<int32_t>(element_, size_);
                    sample_ += delta_;
                }
                else {
                    sample_ = read_bigendian<int32_t>(element_, size_);
                }
                init_ = false;
            }
            return (T)sample_;
        }

        ptrdiff_t operator-(const TraceIterator& rval) const {
            return (element_ - rval.element_) / pitch_;
        }

        bool operator==(const TraceIterator& rval) const {
            return element_ == rval.element_;
        }

        bool operator!=(const TraceIterator& rval) const {
            return element_ != rval.element_;
        }

    private:
        const char* element_;
        size_t size_;
        size_t pitch_;
        int32_t delta_;
        int32_t sample_;
        bool transform_;
        bool init_;
    };

private:
    static const size_t BASE_SIZE_V3 = 9;

#pragma pack(push, 1)

    struct ScfHeader
    {
        uint32_t magic_number;
        BigEndianU32 samples;          /* Number of elements in Samples matrix */
        BigEndianU32 samples_offset;   /* Byte offset from start of file */
        BigEndianU32 bases;            /* Number of bases in Bases matrix */
        BigEndianU32 bases_left_clip;  /* OBSOLETE: No. bases in left clip (vector) */
        BigEndianU32 bases_right_clip; /* OBSOLETE: No. bases in right clip (qual) */
        BigEndianU32 bases_offset;     /* Byte offset from start of file */
        BigEndianU32 comments_size;    /* Number of bytes in Comment section */
        BigEndianU32 comments_offset;  /* Byte offset from start of file */
        char version[4];               /* "version.revision" */
        BigEndianU32 sample_size;      /* Size of samples in bytes 1=8bits, 2=16bits*/
        BigEndianU32 code_set;         /* code set used */
        BigEndianU32 private_size;     /* No. of bytes of Private data, 0 if none */
        BigEndianU32 private_offset;   /* Byte offset from start of file */
        uint32_t spare[18];            /* Unused */

        uint32_t SampleSize() const {
            return version[0] < '2' ? 1 : sample_size;
        }

        bool DeltaTransform() const {
            return version[0] > '2';
        }
    };

    template<typename T>
    struct ScfSamples
    {
        T sample_A;           /* Sample for A trace */
        T sample_C;           /* Sample for C trace */
        T sample_G;           /* Sample for G trace */
        T sample_T;           /* Sample for T trace */
    };

    struct ScfBase
    {
        BigEndianU32 peak_index; /* Index into Samples matrix for base posn */
        uint8_t prob_A;			 /* Probability of it being an A */
        uint8_t prob_C;			 /* Probability of it being an C */
        uint8_t prob_G;			 /* Probability of it being an G */
        uint8_t prob_T;			 /* Probability of it being an T */
        char base;				 /* Called base character */
        uint8_t spare[3];        /* Spare */
    };

    struct ScfPrivate
    {
        uint32_t sig;
        int16_t left;
        int16_t right;
        bool reversed;
    };

#pragma pack(pop)

public:
    // Number of traces in a file.
    static const size_t TRACE_COUNT = 4;

    enum class SearchResult
    {
        SUCCESS,
        NOT_FOUND
    };

    ScfFile(const char* path);

    ScfFile(const ScfFile&) = delete;

    ScfFile(ScfFile&&) = delete;

    ScfFile& operator=(const ScfFile&) = delete;

    ScfFile& operator=(ScfFile&&) = delete;

    bool Sequence(Iterator<char>& begin, Iterator<char>& end) const;

    template <typename T>
    bool Quality(Iterator<T>& begin, Iterator<T>& end) const;

    template <typename T>
    bool Traces(char base, TraceIterator<T>& begin, TraceIterator<T>& end) const {
        if (!header_->samples)
            return false;
        uint32_t index = BaseIndex(base);
        if (index == ~0u)
            return false;
        uint32_t sample_size = header_->SampleSize();
        if (!header_->DeltaTransform()) {
            begin = TraceIterator<T>(file_buffer_.get() + header_->samples_offset + index * sample_size, sample_size, sample_size * TRACE_COUNT, false);
            end = TraceIterator<T>(file_buffer_.get() + header_->samples_offset + header_->samples * sample_size * TRACE_COUNT + index * sample_size, sample_size, sample_size * TRACE_COUNT, false);
        }
        else {
            begin = TraceIterator<T>(file_buffer_.get() + header_->samples_offset + index * header_->samples * sample_size, sample_size, sample_size, true);
            end = TraceIterator<T>(file_buffer_.get() + header_->samples_offset + (index + 1) * header_->samples * sample_size, sample_size, sample_size, true);
        }
        return true;
    }

    template <typename T>
    bool Peaks(Iterator<T>& begin, Iterator<T>& end) const;

    template <typename T>
    SearchResult SearchTag(const char* tag, int32_t number, Iterator<T>& begin, Iterator<T>& end) const;
    
    SearchResult Trims(int32_t& left, int32_t& right) const;

    bool Reversed() const;

private:
    void ValidateHeader();

    void DeltaDecodeTraces();

    static uint32_t BaseIndex(char base) {
        base = LookupTables::Uppercase(base);
        for (uint32_t i = 0; i < TRACE_COUNT; ++i)
            if (base == base_order[i])
                return i;
        return ~0u;
    }

    static const char base_order[TRACE_COUNT];

    std::unique_ptr<char[]> file_buffer_;
    size_t file_size_;
    const ScfHeader* header_;
};