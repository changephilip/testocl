#include "gene.h"
#include <sys/stat.h>
#include "cl_cpp_utility.cpp"
#include "parse.h"
BOOST_COMPUTE_ADAPT_STRUCT(Junction, Junction, (start_, end_))
BOOST_COMPUTE_ADAPT_STRUCT(read_core_t, read_core_t, (junctionCount, junctions))
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
                start_=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfRead, NULL);
                end_=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                sizeof(uint64_t) * numOfRead, NULL);
                strand=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint8_t) * numOfRead, NULL);
                core=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                sizeof(read_core_t) * numOfRead, NULL);
        }
};
cl_d_Reads cl_chipMallocRead(cl::CommandQueue &queue, cl::Context &contexts,
                             h_Reads &h_reads, int32_t numOfRead)
{
        cl_d_Reads cl_d_reads(contexts, numOfRead);
        queue.enqueueWriteBuffer(cl_d_reads.start_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfRead,
                                 &(h_reads.start_[0]));
        queue.enqueueWriteBuffer(cl_d_reads.end_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfRead,
                                 &(h_reads.end_[0]));
        queue.enqueueWriteBuffer(cl_d_reads.strand, CL_FALSE, 0,
                                 sizeof(uint8_t) * numOfRead,
                                 &(h_reads.strand[0]));
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
                 start_=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfBin, NULL);
                 end_=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                sizeof(uint64_t) * numOfBin, NULL);
                 strand=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint8_t) * numOfBin, NULL);
                 core=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                sizeof(bin_core_t) * numOfBin, NULL);
        }
};

cl_d_Bins cl_chipMallocBin(cl::CommandQueue &queue, cl::Context &contexts,
                           h_Bins &h_bins, int32_t numOfBin)
{
        cl_d_Bins cl_d_bins(contexts, numOfBin);
        queue.enqueueWriteBuffer(cl_d_bins.start_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfBin,
                                 &(h_bins.start_[0]));
        queue.enqueueWriteBuffer(cl_d_bins.end_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfBin,
                                 &(h_bins.end_[0]));
        queue.enqueueWriteBuffer(cl_d_bins.strand, CL_FALSE, 0,
                                 sizeof(uint8_t) * numOfBin,
                                 &(h_bins.strand[0]));
        queue.enqueueWriteBuffer(cl_d_bins.core, CL_FALSE, 0,
                                 sizeof(bin_core_t) * numOfBin,
                                 &(h_bins.core[0]));
        return cl_d_bins;
}
struct cl_d_ASEs {
        cl::Buffer start_;
        cl::Buffer end_;
        cl::Buffer strand;

        cl::Buffer core;

        cl_d_ASEs(cl::Context &contexts, int32_t numOfASE)
        {
		cl_int err;
                 start_=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint64_t) * numOfASE,NULL, &err);
		checkCLBuffer(err,__LINE__);
                 end_=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                sizeof(uint64_t) * numOfASE,NULL, &err);
		checkCLBuffer(err,__LINE__);
                 strand=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                  sizeof(uint8_t) * numOfASE,NULL, &err);
		checkCLBuffer(err,__LINE__);
                 core=cl::Buffer(contexts, CL_MEM_READ_WRITE,
                                sizeof(ase_core_t) * numOfASE,NULL, &err);
        }
};

cl_d_ASEs cl_chipMallocASE(cl::CommandQueue &queue, cl::Context &contexts,
                           h_ASEs &h_ases, int32_t numOfASE)
{
        cl_d_ASEs cl_d_ases(contexts, numOfASE);
	cl_int err;
        err = queue.enqueueWriteBuffer(cl_d_ases.start_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfASE,
                                 &(h_ases.start_[0]));
	checkCLIOBuffer(err,__LINE__);
        queue.enqueueWriteBuffer(cl_d_ases.end_, CL_FALSE, 0,
                                 sizeof(uint64_t) * numOfASE,
                                 &(h_ases.end_[0]));
        queue.enqueueWriteBuffer(cl_d_ases.strand, CL_FALSE, 0,
                                 sizeof(uint8_t) * numOfASE,
                                 &(h_ases.strand[0]));
        queue.enqueueWriteBuffer(cl_d_ases.core, CL_FALSE, 0,
                                 sizeof(ase_core_t) * numOfASE,
                                 &(h_ases.core[0]));
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
        int32_t numOfRead = int32_t(h_reads.start_.size());
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
        /*
        boost::compute::vector<read_core_t> bc_d_cores((uint64_t)numOfRead,
                                                       bc_context);
        */

        boost::compute::copy(h_reads.start_.begin(), h_reads.start_.end(),
                             bc_d_starts.begin(), bc_queue);
        boost::compute::copy(h_reads.end_.begin(), h_reads.end_.end(),
                             bc_d_ends.begin(), bc_queue);
        boost::compute::copy(h_reads.strand.begin(), h_reads.strand.end(),
                             bc_d_strand.begin(), bc_queue);
        /*
        boost::compute::copy(h_reads.core.begin(), h_reads.core.end(),
                             bc_d_cores.begin(), bc_queue);
        */
        boost::compute::vector<int> d_indices((uint64_t)numOfRead, bc_context);
        /*
        std::vector<int> h_indices(numOfRead);
        for (int tmp=0;tmp<numOfRead;tmp++){
                h_indices[tmp]=tmp;
        }
        */
        // boost::compute::copy(h_indices.begin(),h_indices.end(),d_indices.begin(),bc_queue);

        boost::compute::iota(d_indices.begin(), d_indices.end(), 0, bc_queue);

        /* DONE CL: thrust sort by key*/
        if (boostSort(bc_d_starts, d_indices, bc_queue)) {
                // boost::compute::gather(d_indices.begin(), d_indices.end(),
                // bc_d_starts, bc_d_starts);
                boost::compute::gather(d_indices.begin(), d_indices.end(),
                                       bc_d_ends.begin(), bc_d_ends.begin(),
                                       bc_queue);
                boost::compute::gather(d_indices.begin(), d_indices.end(),
                                       bc_d_strand.begin(), bc_d_strand.begin(),
                                       bc_queue);
                /*
                boost::compute::gather(d_indices.begin(), d_indices.end(),
                                       bc_d_cores.begin(), bc_d_cores.begin(),
                                       bc_queue);
                */
        } else {
                perror("boostSort return False");
                exit(EXIT_FAILURE);
        }
        // copy back to h_reads
        boost::compute::copy(h_reads.end_.begin(), h_reads.end_.end(),
                             bc_d_ends.begin(), bc_queue);
        boost::compute::copy(h_reads.strand.begin(), h_reads.strand.end(),
                             bc_d_strand.begin(), bc_queue);
        /*
        boost::compute::copy(h_reads.core.begin(), h_reads.core.end(),
                             bc_d_cores.begin(), bc_queue);
        */

        /* Start CL setup*/
        // cl::Contexts contexts = initCppCLSetup();
        env CLEnv = initCppCLSetup();
	cl_int err;
        cl::CommandQueue queue(CLEnv.context, CLEnv.selectedDevices[0], 0,
                               NULL);
        /* import CL kernel first*/
        std::string kernelFileName = "";
        // kernelList allKernel = initCompileKernel_List(
        //   CLEnv.selectedDevices, CLEnv.context, kernelFileName);
        kernelList allKernel =
            initCompileKernel_List(CLEnv.selectedDevices, CLEnv.context);

        cl::Buffer cl_assist_reads(CLEnv.context, CL_MEM_READ_WRITE,
                                   sizeof(Assist) * numOfRead, NULL);
        cl::NDRange offset(0, 0);
        cl::NDRange global_size(nBlock*blockSize, 1);
        cl::NDRange local_size(blockSize, 1);

        // allocate memory on gpu
        cl_d_Bins cl_d_bins =
            cl_chipMallocBin(queue, CLEnv.context, h_bins, numOfBin);
        cl_d_Reads cl_d_reads =
            cl_chipMallocRead(queue, CLEnv.context, h_reads, numOfRead);
        cl_d_ASEs cl_d_ases =
            cl_chipMallocASE(queue, CLEnv.context, h_ases, numOfASE);

        // copy and gather bc_d_cores
        // first transfer indices to host
        std::vector<int> h_indices(numOfRead);
        boost::compute::copy(d_indices.begin(),d_indices.end(),h_indices.begin(),bc_queue);
        int32_t *indices;
        indices = new int32_t[numOfRead];
        for (int32_t i = 0; i < numOfRead; ++i) {
                indices[i] = h_indices[i];
        }
        // allocate cl_d_indices, cores_in on device and h_core on host
        cl::Buffer cl_d_indices(CLEnv.context, CL_MEM_READ_WRITE,
                                sizeof(int32_t) * numOfRead, NULL);
        cl::Buffer cl_d_cores_in(CLEnv.context, CL_MEM_READ_WRITE,
                              sizeof(read_core_t) * numOfRead, NULL);
        read_core_t * h_core;
        h_core = new read_core_t [numOfRead];
        for (int32_t i=0; i < numOfRead; ++i) {
                h_core[i] = h_reads.core[i];
        }

	queue.finish();
        err = queue.enqueueWriteBuffer(cl_d_indices, CL_TRUE, 0, sizeof(int32_t)*numOfRead,indices );
	checkCLIOBuffer(err,__LINE__);
        err = queue.enqueueWriteBuffer(cl_d_cores_in,CL_TRUE,0,sizeof(read_core_t)*numOfRead,h_core);
	checkCLIOBuffer(err,__LINE__);
        //cl::Buffer bc_d_cores(CLEnv.context,CL_MEM_READ_WRITE,sizeof(read_core_t)*numOfRead,NULL);
        //gather cores and tranfer to cl_d_reads.core
        autoSetKernelArgs(allKernel.gather_kernel, cl_d_indices,cl_d_cores_in,cl_d_reads.core,(uint64_t)numOfRead);
        err = queue.enqueueNDRangeKernel(allKernel.gather_kernel, offset, global_size, local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();
        delete[] indices;
        delete[] h_core;
        // assign reads to bins
        std::cout << "starting assign reads..." << std::endl;
        /*DONE CL: assign read kernel*/
        autoSetKernelArgs(allKernel.gpu_assign_read_kernel, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, cl_d_reads.start_, cl_d_reads.end_,
                          cl_d_reads.strand, cl_d_reads.core, numOfRead,
                          cl_assist_reads);
        err = queue.enqueueNDRangeKernel(allKernel.gpu_assign_read_kernel, offset,
                                   global_size, local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();

        int32_t tpmSize = (nextPow2(nBlock) + 1) / 2;
        cl::Buffer d_tempTPM(CLEnv.context, CL_MEM_READ_WRITE,
                             sizeof(float) * numOfBin, NULL);
        cl::Buffer d_tpmCounter(CLEnv.context, CL_MEM_READ_WRITE,
                                sizeof(float) * tpmSize, NULL);
        queue.enqueueFillBuffer(d_tempTPM, 0, 0, numOfBin*sizeof(float));
        queue.enqueueFillBuffer(d_tpmCounter, 0, 0, tpmSize*sizeof(float));

        std::cout << "starting count tpm..." << std::endl;
        /*DONE CL:count tempTPM,reduce singlePass,countTPM*/
        autoSetKernelArgs(allKernel.gpu_count_tempTPM, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, d_tempTPM);
        err = queue.enqueueNDRangeKernel(allKernel.gpu_count_tempTPM, offset,
                                   global_size, local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();
        // autoSetKernelArgs(allKernel.gpu_block_reduce, blockSize, d_tempTPM,
        //                         d_tpmCounter);

        // queue.enqueueNDRangeKernel(allKernel.gpu_block_reduce, offset,
        //                              global_size, local_size);
        // queue.enqueueBarrierWithWaitList();
        boost::compute::vector<float> bc_tempTPM(numOfBin, bc_context);
        float bc_tpmCounter = 0.0f;
        boost::compute::reduce(bc_tempTPM.begin(), bc_tempTPM.end(),
                               &bc_tpmCounter,boost::compute::plus<float>(),bc_queue);
        queue.enqueueWriteBuffer(d_tpmCounter, CL_TRUE, 0, sizeof(float),
                                 &bc_tpmCounter);
        autoSetKernelArgs(allKernel.gpu_count_TPM, cl_d_bins.start_,
                          cl_d_bins.end_, cl_d_bins.strand, cl_d_bins.core,
                          numOfBin, d_tempTPM, d_tpmCounter);
        err = queue.enqueueNDRangeKernel(allKernel.gpu_count_TPM, offset, global_size,
                                   local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();

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
        err = queue.enqueueNDRangeKernel(allKernel.gpu_assign_ASE_kernel, offset,
                                   global_size, local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();

        // auxiliary array
        cl::Buffer d_assist_read_ases(CLEnv.context, CL_MEM_READ_WRITE,
                                      sizeof(Assist) * numOfASE);
        cl::Buffer ACT(CLEnv.context, CL_MEM_READ_WRITE,
                       sizeof(ASECounter) * numOfASE);

        // compute number of thread block
        nBlock = (unsigned(numOfASE) + blockSize - 1) / blockSize;
        cl::NDRange global_size_ASE(nBlock*blockSize, 1);
        // assign reads to ases
        /*TODO CL:gpu assign read ASE*/
        autoSetKernelArgs(allKernel.gpu_assign_read_ASE_kernel,
                          cl_d_ases.start_, cl_d_ases.end_, cl_d_ases.strand,
                          cl_d_ases.core, numOfASE, cl_d_reads.start_,
                          cl_d_reads.end_, cl_d_reads.strand, cl_d_reads.core,
                          numOfRead, d_assist_read_ases, ACT);
        err = queue.enqueueNDRangeKernel(allKernel.gpu_assign_read_ASE_kernel, offset,
                                   global_size_ASE, local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();
        // count psi
        /*DONE CL: count PSI*/
        size_t psiSize = sizeof(ASEPsi) * numOfASE;
        cl::Buffer d_ase_psi(CLEnv.context, CL_MEM_READ_WRITE, psiSize);
        std::cout << "starting count psi..." << std::endl;

        autoSetKernelArgs(allKernel.gpu_count_PSI, cl_d_ases.start_,
                          cl_d_ases.end_, cl_d_ases.strand, cl_d_ases.core,
                          numOfASE, d_ase_psi, ACT);

        err = queue.enqueueNDRangeKernel(allKernel.gpu_count_PSI, offset,
                                   global_size_ASE, local_size);
	checkCLEnqueue(err,__LINE__);
	queue.finish();

        ASEPsi *h_ase_psi = new ASEPsi[numOfASE];

        queue.enqueueReadBuffer(d_ase_psi, CL_TRUE, 0, psiSize, h_ase_psi);
        float tpmCounter = 0;
        queue.enqueueReadBuffer(d_tpmCounter, CL_TRUE, 0, sizeof(float),
                                &tpmCounter);

	queue.finish();
	int32_t countIn, countOut;
        UMAP::const_iterator gid, bin;
        countIn = countOut = 0;
        for (int i = 0; i < numOfASE; i++) {
                if (h_ase_psi[i].countIn)
                        countIn++;
                if (h_ase_psi[i].countOut)
                        countOut++;
                gid = g_gid_map.find(h_ase_psi[i].gid_h);
                bin = g_name_map.find(h_ase_psi[i].bin_h);
                std::cout << gid->second << " "
                          << ((bin == g_name_map.end()) ? "null" : bin->second)
                          << " " << h_ase_psi[i].countIn << " "
                          << h_ase_psi[i].countOut << " " << h_ase_psi[i].psi
                          << std::endl;
        }
        std::cout << countIn << " " << countOut << std::endl;
        // tpm
        std::cout << "tpm count: " << tpmCounter << std::endl;
        // free memory
        std::cout << "free memory..." << std::endl;
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

        cl_HandleBin(h_bins, h_reads, h_ases);

        gettimeofday(&t2_time, 0);
        std::cout << "computing spent time: "
                  << (float)(t2_time.tv_sec - t_time.tv_sec) << "s"
                  << std::endl;
        return 0;
}
