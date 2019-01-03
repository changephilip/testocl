#ifndef _CHIP_BIN_CUH_
#define _CHIP_BIN_CUH_

#define SE_ANCHOR
#ifdef SE_ANCHOR
#define anchorCount 4
#define coordinateCount 6
// serialization
#else
#ifdef RI_ANCHOR
constant int anchorCount = 2;
constant int coordinateCount = 4;
#endif
#endif

//#define CL_ALIGNED(_x) __attribute__((aligned(_x)))
typedef unsigned char uint8_t;
typedef unsigned long long uint64_t;
typedef int int32_t;
typedef unsigned int uint32_t;
// maximum field size
constant int nameSize = 24;
constant int gidSize = 96;

#define SINGLE_END
#ifdef SINGLE_END
#define junctionSize 5
#elif defined(PAIR_END)
constant int junctionSize = 10;
#endif

// constant uint64_t refLength = (uint64_t)2 << 31;
constant uint64_t refLength = 0x100000000;
// constant uint64_t invalidLength = refLength << 5;
constant uint64_t invalidLength = 0x2000000000;

// kernel parameters
constant int blockSize = 1024;

// psi model
constant float step_psi = 0.01;
constant int readLength = 100;

typedef struct CL_ALIGNED(4) {
    // int32_t start_ = 0;
    // int32_t end_ = 0;
    int32_t start_;
    int32_t end_;

} Junction, Anchor, Assist;

typedef struct CL_ALIGNED(8)  {
    // with junction
    // uint32_t junctionCount = 0;
    uint32_t junctionCount;
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
    //        for (int i = 0; i < junctionSize; i++) junctions[i].start_ =
    //        junctions[i].end_ = 0;
    //    }
}read_core_t;

typedef struct {
    int32_t gid_h;
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

typedef struct CL_ALIGNED(32)  {
    Anchor artRange;
    // int32_t anchor[anchorCount] = {0};
    int32_t anchor[anchorCount];

    /*
    __host__ __device__
    ASECounter() {
        c_memset(&anchor, 0, sizeof(int32_t) * anchorCount);
    }
    */
}ASECounter;

typedef struct CL_ALIGNED(32) {
    // int32_t name_h = 0;
    // uint32_t readCount = 0;
    // float tpmCount = 0.0;
    int32_t name_h;
    uint32_t readCount;
    float tpmCount;

    // uint32_t aseCount;

    /*
    __host__ __device__
    bin_core_t(int32_t name) {
        readCount = tpmCount = 0;
        name_h = name;
    }

    bin_core_t() {}
    */
}bin_core_t ;

typedef struct  {
    // int32_t gid_h = 0;
    // int32_t bin_h = 0;
    int32_t gid_h;
    int32_t bin_h;

    int32_t coordinates[coordinateCount];
    /*
__host__ __device__
ase_core_t(int32_t gid) { gid_h = gid; bin_h = 0; }

ase_core_t() {}
    */
}ase_core_t;

typedef struct {
    int32_t gid_h;
    int32_t bin_h;
    float countIn;
    int32_t countOut;
    float psi;  // by sum(anchors)
    // confidence interval
    float ciStart;
    float ciEnd;
} ASEPsi;
ASEPsi iniASEPsi(int32_t gid_h,int32_t bin_h,float countIn,int32_t countOut,float psi,float ciStart,float ciEnd){
	ASEPsi thisASEPsi = {gid_h,bin_h,countIn,countOut,psi,ciStart,ciEnd};
	return thisASEPsi; 
}
__kernel void gather_kernel(__global int32_t *indices,
                            __global read_core_t *d_cores_in,
                            __global read_core_t *d_reads_core,
                            uint64_t numOfRead)
{
    int32_t threadId = get_group_id(0) * get_local_size(0) + get_local_id(0);
    if (threadId < numOfRead) {
        d_reads_core[threadId] = d_cores_in[indices[threadId]];
    }
}

void gpu_try_assign_kernel(uint64_t bin_start, uint64_t bin_end, uint32_t id,
                           __global uint64_t *d_starts, int32_t numOfEntry,
                           __global Assist *d_assist)
{
    uint64_t start;
    int32_t left = 0, right = numOfEntry, mid_l, mid_r;

    while (left < right) {
        mid_l = (left + right) / 2;
        start = d_starts[mid_l];
        if (start < bin_start)
            left = mid_l + 1;
        else
            right = mid_l;
    }
    if (left != numOfEntry)
        d_assist[id].start_ = left;
    else {
        d_assist[id].start_ = d_assist[id].end_ = 0;
        return;
    }
    left = 0;
    right = numOfEntry;
    while (left < right) {
        mid_r = (left + right) / 2;
        start = d_starts[mid_r];
        if (start < bin_end)
            left = mid_r + 1;
        else
            right = mid_r;
    }
    if (left)
        d_assist[id].end_ = left;
    else {
        d_assist[id].start_ = d_assist[id].end_ = 0;
    }
}

__kernel void gpu_assign_read_kernel(
    __global uint64_t *d_bins_start_, __global uint64_t *d_bins_end_,
    __global uint8_t *d_bins_strand, __global bin_core_t *d_bins_core,
    int32_t numOfBin, __global uint64_t *d_reads_start_,
    __global uint64_t *d_reads_end_, __global uint8_t *d_reads_strand,
    __global read_core_t *d_reads_core, int32_t numOfRead,
    __global Assist *d_assist)
{
    int32_t binId = get_group_id(0) * get_local_size(0) + get_local_id(0);
    int temp = 0;
    if (binId < numOfBin) {
        gpu_try_assign_kernel(d_bins_start_[binId], d_bins_end_[binId], binId,
                              d_reads_start_, numOfRead, d_assist);
    };
    mem_fence(CLK_GLOBAL_MEM_FENCE);

    for (int readId = d_assist[binId].start_; readId < d_assist[binId].end_;
         readId++) {
        if ((d_reads_strand[readId] != d_bins_strand[binId]) ||
            (d_reads_end_[readId] > d_bins_end_[binId]))
            temp++;
    }
    d_bins_core[binId].readCount =
        d_assist[binId].end_ - d_assist[binId].start_ - temp;
}

__kernel void gpu_count_tempTPM(__global uint64_t *d_bins_start_,
                                __global uint64_t *d_bins_end_,
                                __global uint8_t *d_bins_strand,
                                __global bin_core_t *d_bins_core,
                                int32_t numOfBin, __global float *d_tempTPM)
{
    int32_t binId = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (binId < numOfBin) {
        d_tempTPM[binId] = (float)(d_bins_core[binId].readCount) /
                           (float)(d_bins_end_[binId] - d_bins_start_[binId]);
    }
}

__kernel void gpu_count_TPM(__global uint64_t *d_bins_start_,
                            __global uint64_t *d_bins_end_,
                            __global uint8_t *d_bins_strand,
                            __global bin_core_t *d_bins_core,
                            int32_t numOfBin, __global float *d_tempTPM,
                            __global float *d_tpmCounter)
{
    int32_t binId = get_group_id(0) * get_local_size(0) + get_local_id(0);
    if (binId < numOfBin) {
        if (*d_tpmCounter == 0) return;
        d_bins_core[binId].tpmCount =
            1000000 * d_tempTPM[binId] / (*d_tpmCounter);
    }
}

__kernel void gpu_assign_ASE_kernel(
    __global uint64_t *d_bins_start_, __global uint64_t *d_bins_end_,
    __global uint8_t *d_bins_strand, __global bin_core_t *d_bins_core,
    int32_t numOfBin, __global uint64_t *d_ases_start_,
    __global uint64_t *d_ases_end_, __global uint8_t *d_ases_strand,
    __global ase_core_t *d_ases_core, int32_t numOfASE,
    __global Assist *d_assist)
{
    int32_t binId = get_group_id(0) * get_local_size(0) + get_local_id(0);
    if (binId < numOfBin) {
        gpu_try_assign_kernel(d_bins_start_[binId], d_bins_end_[binId], binId,
                              d_ases_start_, numOfASE, d_assist);
        mem_fence(CLK_GLOBAL_MEM_FENCE);
    }

    for (int32_t aseId = d_assist[binId].start_; aseId < d_assist[binId].end_;
         aseId++) {
        if ((d_ases_strand[aseId] == d_bins_strand[binId]) &&
            (d_ases_end_[aseId] <= d_bins_end_[binId])) {
            d_ases_core[aseId].bin_h = d_bins_core[binId].name_h;
        } else {
            d_ases_core[aseId].bin_h = 0;
        }
    }
}

__kernel void gpu_assign_read_ASE_kernel(
    __global uint64_t *d_ases_start_, __global uint64_t *d_ases_end_,
    __global uint8_t *d_ases_strand, __global ase_core_t *d_ases_core,
    int32_t numOfASE, __global uint64_t *d_reads_start_,
    __global uint64_t *d_reads_end_, __global uint8_t *d_reads_strand,
    __global read_core_t *d_reads_core, int32_t numOfRead,
    __global Assist *d_assist, __global ASECounter *ACT)
{
    int32_t aseId = get_group_id(0) * get_local_size(0) + get_local_id(0);
    uint32_t read_strand, ase_strand, junctionCount;
    uint32_t read_s, read_e, junction_s, junction_e;
    global int32_t *coord;

    if (aseId < numOfASE) {
        gpu_try_assign_kernel(d_ases_start_[aseId], d_ases_end_[aseId], aseId,
                              d_reads_start_, numOfRead, d_assist);
        mem_fence(CLK_GLOBAL_MEM_FENCE);
        coord = d_ases_core[aseId].coordinates;
        ACT[aseId].artRange.start_ = coord[2];
        ACT[aseId].artRange.end_ = coord[3];

        for (int32_t readId = d_assist[aseId].start_;
             readId < d_assist[aseId].end_; readId++) {
            read_strand = d_reads_strand[readId];
            ase_strand = d_ases_strand[aseId];
            if (read_strand == ase_strand) {
                read_s = (uint32_t)(d_reads_start_[readId] & (refLength - 1));
                read_e = (uint32_t)(d_reads_end_[readId] & (refLength - 1));
#ifdef SE_ANCHOR
                junctionCount = d_reads_core[readId].junctionCount;
                if (junctionCount) {
//#pragma unroll
                    for (int32_t jId = 0; jId < junctionCount; jId++) {
                        junction_s =
                            d_reads_core[readId].junctions[jId].start_ +
                            read_s - 1;
                        junction_e =
                            d_reads_core[readId].junctions[jId].end_ + read_s;
                        if (ase_strand) {
                            if (junction_s == coord[1] &&
                                junction_e == coord[2])
                                ACT[aseId].anchor[0]++;
                            if (junction_s == coord[3] &&
                                junction_e == coord[4])
                                ACT[aseId].anchor[1]++;
                            if (junction_s == coord[1] &&
                                junction_e == coord[4])
                                ACT[aseId].anchor[2]++;
                        } else {
                            if (junction_s == coord[5] &&
                                junction_e == coord[2])
                                ACT[aseId].anchor[0]++;
                            if (junction_s == coord[3] &&
                                junction_e == coord[0])
                                ACT[aseId].anchor[1]++;
                            if (junction_s == coord[5] &&
                                junction_e == coord[0])
                                ACT[aseId].anchor[2]++;
                        }
                    }
                } else {
                    if ((read_s >= coord[2] && read_e <= coord[3]) ||
                        (read_e >= coord[2] && read_e <= coord[3])) {
                        ACT[aseId].anchor[3]++;
                    }
                }
#elif defined(RI_ANCHOR)
                junctionCount = d_reads_core[readId].junctionCount;
                if (junctionCount) {
//#pragma unroll
                    for (int32_t jId = 0; jId < junctionCount; jId++) {
                        junction_s = d_reads_core[readId].junctions[jId].start +
                                     +read_s - 1;
                        junction_e =
                            d_reads_core[readId].junctions[jId].end_ + read_s;
                        if (ase_strand) {
                            if (junction_s == coord[1] &&
                                junction_e == coord[2])
                                ACT[aseId].anchor[0]++;
                        } else {
                            if (junction_s == coord[2] &&
                                junction_e == coord[1])
                                ACT[aseId].anchor[0]++;
                        }
                    }
                } else {
                    if (ase_strand) {
                        if ((read_s >= coord[1] && read_s <= coord[2]) ||
                            (read_e >= coord[1] && read_e <= coord[2])) {
                            ACT[aseId].anchor[1]++;
                        }
                    } else {
                        if ((read_s >= coord[2] && read_s <= coord[1]) ||
                            (read_e >= coord[2] && read_e <= coord[1])) {
                            ACT[aseId].anchor[1]++;
                        }
                    }
                }
#endif
            }
        }
    }
}

__kernel void gpu_count_PSI(__global uint64_t *d_ases_start_,
                            __global uint64_t *d_ases_end_,
                            __global uint8_t *d_ases_strand,
                            __global ase_core_t *d_ases_core,
                            int32_t numOfASE, __global ASEPsi *d_ase_psi,
                            __global ASECounter *ACT)
{
    int32_t aseId = get_group_id(0) * get_local_size(0) + get_local_id(0);
    int32_t countOut;
    float countIn, psi;
    ASECounter act;

    if (aseId < numOfASE) {
        act = ACT[aseId];
#ifdef SE_ANCHOR
        countIn =
            act.anchor[0] + act.anchor[1] +
            act.anchor[3] / (float)(act.artRange.end_ - act.artRange.start_);
        countOut = act.anchor[2];
        if (act.anchor[3]) {
            psi = (countIn / 3) / (countIn / 3 + countOut);
        } else {
            psi = (countIn / 2) / (countIn / 2 + countOut);
        }
#elif defined(RI_ANCHOR)
        countIn = (float)(act.anchor[0]);
        countOut = act.anchor[1];
        psi = countIn / (countIn + countOut);
#endif
        d_ase_psi[aseId] = iniASEPsi(d_ases_core[aseId].gid_h,
                                  d_ases_core[aseId].bin_h,
                                  countIn,
                                  countOut,
                                  psi,
                                  0,
                                  0);
        mem_fence(CLK_GLOBAL_MEM_FENCE);
    }
}
#endif
