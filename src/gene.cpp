#include "gene.h"
#include "cl_prepare.h"
#include "parse.h"
#include <CL/cl.h>
#include <CL/cl.hpp>
#include <boost/compute.hpp>
#include <sys/stat.h>


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
        queue.enqueueWriteBuffer(cl_d_reads.core, CL_FALSE, 0,
                                 sizeof(read_core_t) * numOfRead,
                                 h_reads.core.data());
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


template <typename TKeys, typename TValues>
bool boostSort(boost::compute::vector<TKeys> &keys,
               boost::compute::vector<TValues> &values,
               boost::compute ::command_queue queue)
{
        boost::compute::sort_by_key(keys.begin(), keys.end(), values.begin(),
                                    queue);
        return boost::compute::is_sorted(keys.begin(), keys.end(), queue);
}

void cl_HandleBin(h_Bins &h_bins, h_Reads &h_reads, h_ASEs &h_ases)
{
        /*initialize boost compute env*/
        boost::compute::device bc_gpu =
            boost::compute::system::default_device();
        boost::compute::context bc_context(bc_gpu);
        boost::compute::command_queue bc_queue(bc_context, bc_gpu);

        int32_t numOfBin = int32_t(h_bins.start_.size());
        int32_t numOfRead = int32_t(h_bins.start_.size());
        int32_t numOfASE = int32_t(h_ases.start_.size());
        std::cout << "numOfBin: " << numOfBin << std::endl;
        std::cout << "numOfRead: " << numOfRead << std::endl;
        std::cout << "numOfASE: " << numOfASE << std::endl;
        // compute number of thread block
        unsigned nBlock = (unsigned(numOfBin) + blockSize - 1) / blockSize;

        /*
  member of d_reads only can be used as type cl_mem or cl::memory,so it can't be
  accesible by boost:compute. We malloc new boost::compute::vector to sort and
  gather,and then reuse cl_chipMalloc to memcpy.
 */
        /* DONE CL:allocate boost compute vectors in bc_context for
         * device_starts,ends,strands and cores, and then copy memory*/
        std::cout << "starting sorting reads..." << std::endl;
        boost::compute::vector<uint64_t> bc_d_starts((uint64_t)numOfRead,
                                                     bc_context);
        boost::compute::vector<uint64_t> bc_d_ends((uint64_t)numOfRead,
                                                   bc_context);
        boost::compute::vector<uint8_t> bc_d_strand((uint64_t)numOfRead,
                                                    bc_context);
        boost::compute::vector<read_core_t> bc_d_cores((uint64_t)numOfRead,
                                                       bc_context);

        boost::compute::copy(h_bins.start_.begin(), h_bins.start_.end(),
                             bc_d_starts.begin(), bc_queue);

        boost::compute::vector<int> d_indices((uint64_t)numOfRead, context);
        boost::compute::iota(d_indices.begin(), d_indices.end(), 0, bc_queue);

        /* DONE CL: thrust sort by key*/
        if (boostSort(bc_d_starts, d_indices, bc_queue)) {
                // boost::compute::gather(d_indices.begin(), d_indices.end(),
                // bc_d_starts, bc_d_starts);
                boost::compute::gather(d_indices.begin(), d_indices.end(),
                                       bc_d_ends, bc_d_ends);
                boost::compute::gather(d_indices.begin(), d_indices.end(),
                                       bc_d_strand, bc_d_strand);
                boost::compute::gather(d_indices.begin(), d_indices.end(),
                                       bc_d_cores, bc_d_cores);
        } else {
                perror("boostSort return False");
                exit(EXIT_FAILURE);
        }
        // copy back to h_bins
        boost::compute::copy(h_bins.end_.begin(), h_bins.end_.end(),
                             bc_d_ends.begin(), bc_queue);
        boost::compute::copy(h_bins.strand.begin(), h_bins.strand.end(),
                             bc_d_strand.begin(), bc_queue);
        boost::compute::copy(h_bins.core.begin(), h_bins.core.end(),
                             bc_d_cores.begin(), bc_queue);

        /* Start CL setup*/
        // cl::Contexts contexts = initCppCLSetup();
        env CLEnv = initCppCLSetup();
        cl::CommandQueue queue(CLEnv.context, CLEnv.selectedDevices[0], 0,
                               NULL);
        /* import CL kernel first*/
        std::string kernelFileName = "";
        kernelList allKernel = initCompileKernel_List(
            CLEnv.selectedDevices, CLEnv.context, kernelFileName);
        cl::Buffer cl_assist_reads(CLEnv.context, CL_MEM_READ_WRITE,
                                   sizeof(Assist) * numOfRead, NULL);
        cl::NDRange offset(0, 0);
        cl::NDRange global_size(numOfBin, 1);
        cl::NDRange local_size(blockSize, 1);

        // allocate memory on gpu
        cl_d_Bins cl_d_bins =
            cl_chipMallocBin(queue, CLEnv.context, h_bins, numOfBin);
        cl_d_Reads cl_d_reads =
            cl_chipMallocRead(queue, CLEnv.context, h_reads, numOfRead);
        cl_d_ASEs cl_d_ases =
            cl_chipMallocASE(queue, CLEnv.context, h_ases, numOfASE);

        // assign reads to bins
        std::cout << "starting assign reads..." << std::endl;
        /*TODO CL: assign read kernel*/
        autoSetKernelArgs(allKernel.gpu_assign_read_kernel, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, cl_d_reads.start_, cl_d_reads.end_,
                          cl_d_reads.strand, cl_d_reads.core, numOfRead,
                          cl_assist_reads);
        queue.enqueueNDRangeKernel(allKernel.gpu_assign_read_kernel, offset,
                                   global_size, local_size);
        queue.enqueueBarrierWithWaitList();

        int32_t tpmSize = (nextPow2(nBlock) + 1) / 2;
        cl::Buffer d_tempTPM(CLEnv.context, CL_MEM_READ_WRITE,
                             sizeof(float) * numOfBin, NULL);
        cl::Buffer d_tpmCounter(CLEnv.context, CL_MEM_READ_WRITE,
                                sizeof(float) * tpmSize, NULL);
        queue.enqueueFillBuffer(d_tempTPM, 0, 0, numOfBin);
        queue.enqueueFillBuffer(d_tpmCounter, 0, 0, tpmSize);

        std::cout << "starting count tpm..." << std::endl;
        /*TODO CL:count tempTPM,reduce singlePass,countTPM*/
        autoSetKernelArgs(allKernel.gpu_count_tempTPM, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, d_tempTPM);
        queue.enqueueNDRangeKernel(allKernel.gpu_count_tempTPM, offset,
                                   global_size, local_size);
        queue.enqueueBarrierWithWaitList();
        // autoSetKernelArgs(allKernel.gpu_block_reduce, blockSize, d_tempTPM,
        //                         d_tpmCounter);

        // queue.enqueueNDRangeKernel(allKernel.gpu_block_reduce, offset,
        //                              global_size, local_size);
        //queue.enqueueBarrierWithWaitList();
        boost::compute::vector<float> bc_tempTPM(numOfBin,bc_context);
        float bc_tpmCounter=0.0f;
        boost::compute::reduce(bc_tempTPM.begin(),bc_tempTPM.end(),&bc_tpmCounter);
        queue.enqueueWriteBuffer(d_tpmCounter, CL_TRUE, sizeof(Float),d_tpmCounter);
        autoSetKernelArgs(allKernel.gpu_count_TPM, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, d_tempTPM, d_tpmCounter);
        queue.enqueueNDRangeKernel(allKernel.gpu_count_TPM, offset, global_size,
                                   local_size);
        queue.enqueueBarrierWithWaitList();

        // auxiliary array
        cl::Buffer d_assist_ases(CLEnv.context, CL_MEM_READ_WRITE,
                                 sizeof(Assist) * numOfBin);
        // assign ases to bins
        /*DONE CL:assign ASE kernel*/
        std::cout << "starting assign ases..." << std::endl;
        autoSetKernelArgs(allKernel.gpu_assign_ASE_kernel, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, cl_d_ases.start_, cl_d_ases.end_,
                          cl_d_ases.strand, cl_d_ases.core, numOfASE,
                          d_assist_ases);
        queue.enqueueNDRangeKernel(allKernel.gpu_assign_ASE_kernel, offset,
                                   global_size, local_size);
        queue.enqueueBarrierWithWaitList();

        // auxiliary array
        cl::Buffer d_assist_read_ases(CLEnv.context, CL_MEM_READ_WRITE,
                                      sizeof(Assist) * numOfASE);
        cl::Buffer ACT(CLEnv.Context, CL_MEM_READ_WRITE,
                       sizeof(ASECounter) * numOfASE);

        // compute number of thread block
        nBlock = (unsigned(numOfASE) + blockSize - 1) / blockSize;
        cl::NDRange global_size_ASE(nBlock, 1);
        // assign reads to ases
        /*TODO CL:gpu assign read ASE*/
        autoSetKernelArgs(allKernel.gpu_assign_read_ASE_kernel,
                          cl_d_ases.start_, cl_d_ases.end_, cl_d_ases.strand,
                          cl_d_ases.core, numOfASE, cl_d_reads.start_,
                          cl_d_reads.end_, cl_d_reads.strand, cl_d_reads.core,
                          numOfRead, d_assist_read_ases, ACT);
        queue.enqueueNDRangeKernel(allKernel.gpu_assign_read_ASE_kernel, offset,
                                   global_size_ASE, local_size);
        queue.enqueueBarrierWithWaitList();
        // count psi
        /*TODO CL: count PSI*/
        size_t psiSize = sizeof(ASEPsi) * numOfASE;
        cl::Buffer d_ase_psi(CLEnv.context, CL_MEM_READ_WRITE, psiSize);
        std::cout << "starting count psi..." << std::endl;

        autoSetKernelArgs(allKernel.gpu_count_PSI, cl_d_ases.start_,
                          cl_d_ases.end_, cl_d_ases.strand, cl_d_ases.core,
                          numOfASE, d_ase_psi, ACT);

        queue.enqueueNDRangeKernel(allKernel.gpu_count_PSI, offset,
                                   global_size_ASE, local_size);
        queue.enqueueBarrierWithWaitList();

        ASEPsi *h_ase_psi = new ASEPsi[numOfASE];

        queue.enqueueReadBuffer(d_ase_psi, CL_TRUE, 0, psiSize, h_ase_psi);
        float tpmCounter = 0;
        queue.enqueueReadBuffer(d_tpmCounter, CL_TRUE, 0, sizeof(float),
                                &tpmCounter);

        std::cout << "free memory..." << std::endl;

        // free memory
        delete[] h_ase_psi;
}

int main(int argc, char **argv)
{
        if (argc != 4) {
                std::cerr << "usage: " << argv[0]
                          << " [GFF file path] [BAM file path] [GFF file path]"
                          << std::endl;
                return 1;
        }

        struct timeval start_time, t_time, t2_time;
        gettimeofday(&start_time, 0);

        // load bins and ases.
        h_Bins h_bins;
        h_ASEs h_ases;

        struct stat buffer;
        if (stat(serFilename, &buffer) == 0) {
                std::cout << "loading bins and ases from serialization..."
                          << std::endl;
                LoadDataFromSerialization(h_bins, h_ases);
        } else {
                std::cout << "loading bins..." << std::endl;
                LoadBinFromGff(h_bins, argv[1]);

                std::cout << "loading ases..." << std::endl;
                LoadAseFromGFF(h_ases, argv[3]);

                std::cout << "saving bins and ases to serialization..."
                          << std::endl;
                SaveDataToSerialization(h_bins, h_ases);
        }

        // load reads
        std::cout << "loading reads..." << std::endl;
        h_Reads h_reads;
        LoadReadFromBam(h_reads, argv[2]);

        gettimeofday(&t_time, 0);
        std::cout << "load spent time: "
                  << (float)(t_time.tv_sec - start_time.tv_sec) << "s"
                  << std::endl;

        std::cout << "start kernel program..." << std::endl;

        cl_HandleBin(h_bins, h_bins, h_ases);

        gettimeofday(&t2_time, 0);
        std::cout << "computing spent time: "
                  << (float)(t2_time.tv_sec - t_time.tv_sec) << "s"
                  << std::endl;
        return 0;
}

