#ifndef _CHIP_HEADER_H_
#define _CHIP_HEADER_H_

#define CL_TARGET_OPENCL_VERSION 120
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <CL/cl.h>
#include <CL/cl.hpp>
#include "c_string.h"
// boost
#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/compute.hpp>
#define _LINUX_
#ifdef _LINUX_
#include <ctime>
#include <sys/time.h>
#endif

#define CL_ALIGN(_x) __attribute__((aligned(_x)))

#define SE_ANCHOR
#ifdef SE_ANCHOR
const int anchorCount = 4;
const int coordinateCount = 6;
// serialization
const char *serFilename = "../gencode_SE.ser";
#else
#ifdef RI_ANCHOR
const int anchorCount = 2;
const int coordinateCount = 4;
// serialization
const char *serFilename = "../gencode_RI.ser";
#endif
#endif

typedef unsigned char uint8_t;

// maximum field size
const int nameSize = 24;
const int gidSize = 96;

#define SINGLE_END
#ifdef SINGLE_END
#define junctionSize  5
#elif defined(PAIR_END)
#define junctionSize  10
#endif

const uint64_t refLength = (uint64_t)2 << 31;
const uint64_t invalidLength = refLength << 5;

// kernel parameters
const int blockSize = 1024;

// psi model
const float step = 0.01;
const int readLength = 100;

typedef struct CL_ALIGN(4) {
    int32_t start_;
    int32_t end_;
} Junction, Anchor, Assist;

struct CL_ALIGN(8) read_core_t {
    // with junction
    uint32_t junctionCount;
    Junction junctions[junctionSize];

    read_core_t() {
        junctionCount = 0;
    }

//    __host__ __device__
//    void clear() {
//        junctionCount = 0;
//        for (int i = 0; i < junctionSize; i++) junctions[i].start_ = junctions[i].end_ = 0;
//    }
};

struct h_Reads {
    std::vector<uint64_t > start_;
    std::vector<uint64_t > end_;
    std::vector<uint8_t > strand;
    std::vector<read_core_t > core;
};

struct d_Reads {
    uint64_t *start_;
    uint64_t *end_;
    uint8_t *strand;

    read_core_t *core;
};

typedef struct {
    size_t gid_h;
    uint32_t anchorIndex;
} JunctionTag;

typedef struct {
    Junction range;
    JunctionTag tag;
} JunctionTagASE;

typedef struct {
    Anchor range;
    JunctionTag tag;
} ASERelated;

struct  CL_ALIGN(32) ASECounter {
    Anchor artRange;
    int32_t anchor[anchorCount];

    ASECounter() {
        c_memset(&anchor, 0, sizeof(int32_t) * anchorCount);
    }
};

struct CL_ALIGN(32)  bin_core_t {
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & name_h;
        ar & readCount;
        ar & tpmCount;
    }

    size_t name_h;
    uint32_t readCount;
    float tpmCount;
    // uint32_t aseCount;

    bin_core_t(size_t name) {
        readCount = tpmCount = 0;
        name_h = name;
    }

    bin_core_t() {}
};

struct h_Bins {
    friend class boost::serialization::access;
    template<class Archive>
        void serialize(Archive &ar, const unsigned int version)
    {
        ar & start_;
        ar & end_;
        ar & strand;
        ar & core;
    }

    std::vector<uint64_t > start_;
    std::vector<uint64_t > end_;
    std::vector<uint8_t > strand;
    std::vector<bin_core_t > core;
};

struct d_Bins {
    uint64_t *start_;                          // absolute starting coordinate of bin
    uint64_t *end_;
    uint8_t *strand;

    bin_core_t *core;
};

struct ase_core_t {
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & gid_h;
        ar & bin_h;
        ar & coordinates;
    }

    size_t gid_h;
    size_t bin_h;
    int32_t coordinates[coordinateCount];

    ase_core_t(size_t gid) { gid_h = gid; bin_h = 0; }

    ase_core_t() {}
};

struct h_ASEs {
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & start_;
        ar & end_;
        ar & strand;
        ar & core;
    }

    std::vector<uint64_t > start_;
    std::vector<uint64_t > end_;
    std::vector<uint8_t > strand;
    std::vector<ase_core_t > core;
};

struct d_ASEs {
    uint64_t *start_;
    uint64_t *end_;
    uint8_t *strand;

    ase_core_t *core;
};

typedef struct {
    size_t gid_h;
    size_t bin_h;
    float countIn;
    int32_t countOut;
    float psi;               // by sum(anchors)
    // confidence interval
    float ciStart;
    float ciEnd;
} ASEPsi;


#endif
