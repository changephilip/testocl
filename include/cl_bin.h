//#ifndef _CHIP_BIN_CUH_
//#define _CHIP_BIN_CUH_

//#include <CL/cl.h>
//#include <CL/cl.hpp>


//#define _LINUX_
//#ifdef _LINUX_
//#endif

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

typedef unsigned char cl_uchar;

// maximum field size
const int nameSize = 24;
const int gidSize = 96;

#define SINGLE_END
#ifdef SINGLE_END
const int junctionSize = 5;
#elif defined(PAIR_END)
const int junctionSize = 10;
#endif

const cl_ulong refLength = (cl_ulong)2 << 31;
const cl_ulong invalidLength = refLength << 5;

// kernel parameters
const int blockSize = 1024;

// psi model
const float step = 0.01;
const int readLength = 100;

typedef struct CL_ALIGNED(4) {
    cl_int start_=0;
    cl_int end_=0;
}Junction, Anchor, Assist;

struct CL_ALIGNED(8) read_core_t {
    // with junction
    cl_uint junctionCount;
    Junction junctions[junctionSize];

    /*
    __host__ __device__
    read_core_t() {
        junctionCount = 0;
    }
    */
//    __host__ __device__
//    void clear() {
//        junctionCount = 0;
//        for (int i = 0; i < junctionSize; i++) junctions[i].start_ = junctions[i].end_ = 0;
//    }
};


typedef struct {
    size_t gid_h;
    cl_uint anchorIndex;
} JunctionTag;

typedef struct {
    Junction range;
    JunctionTag tag;
} JunctionTagASE;

typedef struct {
    Anchor range;
    JunctionTag tag;
} ASERelated;

struct CL_ALIGNED(32) ASECounter {
    Anchor artRange;
    cl_int anchor[anchorCount]= {0};

    /*
    __host__ __device__
    ASECounter() {
        c_memset(&anchor, 0, sizeof(cl_int) * anchorCount);
    }
    */
};

struct CL_ALIGNED(32) bin_core_t {
    size_t name_h;
    cl_uint readCount=0;
    float tpmCount=0.0;
    // cl_uint aseCount;

    /*
    __host__ __device__
    bin_core_t(size_t name) {
        readCount = tpmCount = 0;
        name_h = name;
    }
    
    bin_core_t() {}
    */
};


struct ase_core_t {
    size_t gid_h;
    size_t bin_h=0;
    cl_int coordinates[coordinateCount];

        //ase_core_t(size_t gid) { gid_h = gid; bin_h = 0; }

        //ase_core_t() {}
};



typedef struct {
    size_t gid_h;
    size_t bin_h;
    float countIn;
    cl_int countOut;
    float psi;               // by sum(anchors)
    // confidence interval
    float ciStart;
    float ciEnd;
} ASEPsi;




#endif
