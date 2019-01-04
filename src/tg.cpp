#include "gene.h"
#include <sys/stat.h>
#include "cl_cpp_utility.cpp"
#include "parse.h"
int32_t nextPow2(int32_t x)
{
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return ++x;
}
/*simple allocate memory and return device pointer*/
struct cl_d_Reads {
        cl::Buffer start_;
        cl::Buffer end_;
        cl::Buffer strand;

        cl::Buffer core;
        cl_d_Reads(cl::Context &contexts, int32_t numOfRead)
        {
                cl::Buffer start_(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfRead, NULL);
                cl::Buffer end_(contexts, CL_MEM_READ_WRITE,
                                sizeof(uint64_t) * numOfRead, NULL);
                cl::Buffer strand(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfRead, NULL);
                cl::Buffer core(contexts, CL_MEM_READ_WRITE,
                                sizeof(read_core_t) * numOfRead, NULL);
        }
};
cl_d_Reads cl_chipMallocRead(cl::CommandQueue &queue, cl::Context &contexts,
                             h_Reads &h_reads, int32_t numOfRead)
{
        cl_d_Reads cl_d_reads(contexts, numOfRead);
        queue.enqueueWriteBuffer(cl_d_reads.start_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfRead,
                                 h_reads.start_.data());
        queue.enqueueWriteBuffer(cl_d_reads.end_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfRead,
                                 h_reads.end_.data());
        queue.enqueueWriteBuffer(cl_d_reads.strand, CL_FALSE, 0,
                                 sizeof(uint8_t) * numOfRead,
                                 h_reads.strand.data());
        /*
        queue.enqueueWriteBuffer(cl_d_reads.core, CL_FALSE, 0,
                                 sizeof(read_core_t) * numOfRead,
                                 h_reads.core.data());
        */
        return cl_d_reads;
}

struct cl_d_Bins {
        cl::Buffer start_;
        cl::Buffer end_;
        cl::Buffer strand;

        cl::Buffer core;

        cl_d_Bins(cl::Context &contexts, int32_t numOfBin)
        {
                cl::Buffer start_(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfBin, NULL);
                cl::Buffer end_(contexts, CL_MEM_READ_WRITE,
                                sizeof(uint64_t) * numOfBin, NULL);
                cl::Buffer strand(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfBin, NULL);
                cl::Buffer core(contexts, CL_MEM_READ_WRITE,
                                sizeof(bin_core_t) * numOfBin, NULL);
        }
};

cl_d_Bins cl_chipMallocBin(cl::CommandQueue &queue, cl::Context &contexts,
                           h_Bins &h_bins, int32_t numOfBin)
{
        cl_d_Bins cl_d_bins(contexts, numOfBin);
        queue.enqueueWriteBuffer(cl_d_bins.start_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfBin,
                                 h_bins.start_.data());
        queue.enqueueWriteBuffer(cl_d_bins.end_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfBin,
                                 h_bins.end_.data());
        queue.enqueueWriteBuffer(cl_d_bins.strand, CL_FALSE, 0,
                                 sizeof(uint8_t) * numOfBin,
                                 h_bins.strand.data());
        queue.enqueueWriteBuffer(cl_d_bins.core, CL_FALSE, 0,
                                 sizeof(bin_core_t) * numOfBin,
                                 h_bins.core.data());
        return cl_d_bins;
}
struct cl_d_ASEs {
        cl::Buffer start_;
        cl::Buffer end_;
        cl::Buffer strand;

        cl::Buffer core;

        cl_d_ASEs(cl::Context &contexts, int32_t numOfASE)
        {
                cl::Buffer start_(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfASE, NULL);
                cl::Buffer end_(contexts, CL_MEM_READ_WRITE,
                                sizeof(uint64_t) * numOfASE, NULL);
                cl::Buffer strand(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfASE, NULL);
                cl::Buffer core(contexts, CL_MEM_READ_WRITE,
                                sizeof(ase_core_t) * numOfASE, NULL);
        }
};

cl_d_ASEs cl_chipMallocASE(cl::CommandQueue &queue, cl::Context &contexts,
                           h_ASEs &h_ases, int32_t numOfASE)
{
        cl_d_ASEs cl_d_ases(contexts, numOfASE);
        queue.enqueueWriteBuffer(cl_d_ases.start_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfASE,
                                 h_ases.start_.data());
        queue.enqueueWriteBuffer(cl_d_ases.end_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfASE,
                                 h_ases.end_.data());
        queue.enqueueWriteBuffer(cl_d_ases.strand, CL_FALSE, 0,
                                 sizeof(uint8_t) * numOfASE,
                                 h_ases.strand.data());
        queue.enqueueWriteBuffer(cl_d_ases.core, CL_FALSE, 0,
                                 sizeof(ase_core_t) * numOfASE,
                                 h_ases.core.data());
        return cl_d_ases;
}

void cl_HandleBin(const char* kernelFileName)
{
        /*initialize boost compute env*/
                /* Start CL setup*/
        // cl::Contexts contexts = initCppCLSetup();
        env CLEnv = initCppCLSetup();
        cl::CommandQueue queue(CLEnv.context, CLEnv.selectedDevices[0], 0,
                               NULL);
        /* import CL kernel first*/
        //std::string kernelFileName = "";
        // kernelList allKernel = initCompileKernel_List(
        //   CLEnv.selectedDevices, CLEnv.context, kernelFileName);
        kernelList allKernel =
            initCompileKernel_ListG(CLEnv.selectedDevices, CLEnv.context,kernelFileName);

        }

int main(int argc, char **argv)
{
        if (argc != 4) {
                std::cerr << "usage: " << argv[0]
                          << " [GFF file path] [BAM file path] [GFF file path]"
                          << std::endl;
        }
	char* kernelFileName = argv[1];
        std::cout << "start kernel program..." << std::endl;

        cl_HandleBin(kernelFileName);

        return 0;
}
