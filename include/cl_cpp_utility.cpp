#ifndef CL_VERSION_1_2
#define CL_VERSION_1_2
#endif
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#define __CL_ENABLE_EXCEPTIONS
//#include <CL/cl.h>
#include <CL/cl.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

/*

void displayPlatformInfo(cl_platform_id id, cl_platform_info param_name,
                         const char *paramNameAsStr)
{
        cl_int error = 0;
        size_t paramSize = 0;
        error = clGetPlatformInfo(id, param_name, 0, NULL, &paramSize);
        char *moreInfo = (char *)alloca(sizeof(char) * paramSize);
        error = clGetPlatformInfo(id, param_name, paramSize, moreInfo, NULL);

        if (error != CL_SUCCESS) {
                perror("Unable to find any OpenCl platform information");
                return;
        }
        printf("%s: %s\n", paramNameAsStr, moreInfo);
}
void displayDeviceDetails(cl_device_id id, cl_device_info param_name,
                          const char *paramNameAsStr)
{
        cl_int error = 0;
        size_t paramSize = 0;

        error = clGetDeviceInfo(id, param_name, 0, NULL, &paramSize);
        if (error != CL_SUCCESS) {
                perror("Unable to obtain device info for param\n");
        }

        switch (param_name) {
        case CL_DEVICE_TYPE: {
                cl_device_type *devType = (cl_device_type *)alloca(
                    sizeof(cl_device_type) * paramSize);
                error =
                    clGetDeviceInfo(id, param_name, paramSize, devType, NULL);

                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device info for param\n");
                        return;
                }

                switch (*devType) {
                case CL_DEVICE_TYPE_CPU:
                        printf("CPU detected\n");
                        break;
                case CL_DEVICE_TYPE_GPU:
                        printf("GPU detected\n");
                        break;
                case CL_DEVICE_TYPE_ACCELERATOR:
                        printf("Accelerator detected\n");
                        break;
                case CL_DEVICE_TYPE_DEFAULT:
                        printf("default detected\n");
                        break;
                }
        } break;
        case CL_DEVICE_VENDOR_ID:
        case CL_DEVICE_MAX_COMPUTE_UNITS:
        case CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: {
                cl_uint *ret = (cl_uint *)alloca(sizeof(cl_uint) * paramSize);
                error = clGetDeviceInfo(id, param_name, paramSize, ret, NULL);
                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device info for param\n");
                        return;
                }
                switch (param_name) {
                case CL_DEVICE_VENDOR_ID:
                        printf("\tVENDOR ID: 0x%x\n", *ret);
                        break;
                case CL_DEVICE_MAX_COMPUTE_UNITS:
                        printf(
                            "\tMaximum number of parallel compute units: %d\n",
                            *ret);
                        break;
                case CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:
                        printf("\tMaximum dimensions of global/local work-item "
                               "IDs: %d\n",
                               *ret);
                        break;
                }
        } break;
        case CL_DEVICE_MAX_WORK_ITEM_SIZES: {
                cl_uint maxWIDimensions;
                size_t *ret = (size_t *)alloca(sizeof(size_t) * paramSize);
                error = clGetDeviceInfo(id, param_name, paramSize, ret, NULL);
                error =
                    clGetDeviceInfo(id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
                                    sizeof(cl_uint), &maxWIDimensions, NULL);
                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device info for param\n");
                        return;
                }
                printf("\tMaximum number of work-items in each dimension:(");
                for (cl_int i = 0; i < maxWIDimensions; ++i) {
                        printf("%d ", ret[i]);
                }
                printf(" )\n");
        } break;
        case CL_DEVICE_MAX_WORK_GROUP_SIZE: {
                size_t *ret = (size_t *)alloca(sizeof(size_t) * paramSize);
                error = clGetDeviceInfo(id, param_name, paramSize, ret, NULL);
                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device info for param\n");
                        return;
                }
                printf("\tMaximum number of work-items in a work-group: %d\n",
                       *ret);
        } break;
        case CL_DEVICE_NAME:
        case CL_DEVICE_VENDOR: {
                char data[48];
                error = clGetDeviceInfo(id, param_name, paramSize, data, NULL);
                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device name/vendor info for "
                               "param\n");
                        return;
                }
                switch (param_name) {
                case CL_DEVICE_NAME:
                        printf("\tDevice name %s\n", data);
                        break;
                case CL_DEVICE_VENDOR:
                        printf("\tDevice vendor is %s\n", data);
                        break;
                }
        } break;
        case CL_DEVICE_GLOBAL_MEM_CACHE_SIZE: {
                cl_uint *size = (cl_uint *)alloca(sizeof(cl_uint) * paramSize);
                error = clGetDeviceInfo(id, param_name, paramSize, size, NULL);
                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device name/vendor info for "
                               "param\n");
                        return;
                }
                printf("\tDevice global cacheline size: %d bytes\n", (*size));
                break;
        } break;
        case CL_DEVICE_GLOBAL_MEM_SIZE:
        case CL_DEVICE_MAX_MEM_ALLOC_SIZE: {
                cl_ulong *size =
                    (cl_ulong *)alloca(sizeof(cl_ulong) * paramSize);
                error = clGetDeviceInfo(id, param_name, paramSize, size, NULL);
                if (error != CL_SUCCESS) {
                        perror("Unable to obtain device name/vendor info for "
                               "param\n");
                        return;
                }
                switch (param_name) {
                case CL_DEVICE_GLOBAL_MEM_SIZE:
                        printf("\tDevice global mem: %ld mega-bytes\n",
                               (*size) >> 20);
                        break;
                case CL_DEVICE_MAX_MEM_ALLOC_SIZE:
                        printf(
                            "\tDevice max memory allocation: %ld mega-bytes\n",
                            (*size) >> 20);
                        break;
                }
        } break;
        }
}

void displayDeviceInfo(cl_platform_id id, cl_device_type dev_type)
{
        cl_int error = 0;
        cl_uint numOfDevices = 0;

        error = clGetDeviceIDs(id, dev_type, 0, NULL, &numOfDevices);

        if (error != CL_SUCCESS) {
                perror("Unable to obtain any OpenCL compliant device info");
                exit(1);
        }
        cl_device_id *devices =
            (cl_device_id *)alloca(sizeof(cl_device_id) * numOfDevices);

        error = clGetDeviceIDs(id, dev_type, numOfDevices, devices, NULL);

        if (error != CL_SUCCESS) {
                perror("Unable to obtain any OpenCL compliant device info");
                exit(1);
        }

        printf("Number of detected OpenCL devices: %d\n", numOfDevices);

        for (int i = 0; i < numOfDevices; ++i) {
                displayDeviceDetails()
        }
}
cl_context initCLSetup()
{
        cl_platform_id *platforms;
        cl_uint numOfPlatforms;
        cl_int error;

        error = clGetPlatformIDs(0, NULL, &numOfPlatforms);
        if (error < 0) {
                perror("Unable to find any OpenCL platforms");
                exit(1);
        }

        platforms =
            (cl_platform_id *)alloca(sizeof(cl_platform_id) * numOfPlatforms);
        printf("Number of OpenCL Platforms found: %d\n", numOfPlatforms);

        for (cl_uint i = 0; i < numOfPlatforms; ++i) {
                displayPlatformInfo(platforms[i], CL_PLATFORM_PROFILE,
                                    "CL_PLATFORM_PROFILE");
                displayPlatformInfo(platforms[i], CL_PLATFORM_VERSION,
                                    "CL_PLATFORM_VERSION");
                displayPlatformInfo(platforms[i], CL_PLATFORM_NAME,
                                    "CL_PLATFORM_NAME");
                displayPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
                                    "CL_PLATFORM_VENDOR");
                displayPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS,
                                    "CL_PLATFORM_EXTENSIONS");
        }

        //
}

*/


/*template function to generate automatically set kernel args*/
void autoSetKernelArgs_(cl::Kernel& kernel,int &order){}
template <typename Arg1,class... Args>
void autoSetKernelArgs_(cl::Kernel& kernel, int& order,const Arg1& arg1,const Args& ... args)
{
        //int32_t count = sizeof...(args);
	
        //for (int i = 0; i < count; ++i) {
        //        kernel.setArg(i, std::forward<Args>(args)...);
        //}
        cl_int err;
        printf("%d\t",order);
	err=kernel.setArg(order,arg1);
	if (err !=CL_SUCCESS){
		printf("arg %dth is not right!\n",order);
		if (err == CL_INVALID_KERNEL){
			printf("Invalid kernel\n");
		}
		if (err == CL_INVALID_ARG_INDEX){
			printf("Invalid arg index\n");
		}
		if (err == CL_INVALID_ARG_VALUE){
			printf("Invalid arg value\n");
		}	
	}
  	order = order+1;
	autoSetKernelArgs_(kernel,order,args...);
}
template <class... Args>
void autoSetKernelArgs(cl::Kernel &kernel,const Args& ... args){
	int order =0;
	autoSetKernelArgs_(kernel,order,args...);
	printf("\n");
}
struct env {
        cl::Context context;
        std::vector<cl::Device> selectedDevices;
};

/*select all the NVIDIA devices and return contexts*/
env initCppCLSetup()
{
        std::vector<cl::Platform> platforms;
        std::vector<cl::Device> selectedDevices;
        std::string selectedPlatformStringNV("NVIDIA Corporation");
        std::string selectedPlatformStringAMD("Advanced Micro Devices,Inc");
        env returnEnv;
        unsigned long selectedPlatformID = 0;
        {
                cl::Platform::get(&platforms);
                std::vector<std::vector<cl::Device>> allDevices(
                    platforms.size());

                for (unsigned long i = 0; i < platforms.size(); ++i) {
                        std::string hardware;
                        hardware = platforms[i].getInfo<CL_PLATFORM_NAME>();
                        std::cout << "The " << i << "th platform is "
                                  << hardware << std::endl;
                        if (hardware == selectedPlatformStringNV) {
                                selectedPlatformID = i;
                        }
                }
                platforms[selectedPlatformID].getDevices(CL_DEVICE_TYPE_GPU,
                                                         &selectedDevices);
                cl::Context contexts(selectedDevices);
                returnEnv.context = contexts;
                returnEnv.selectedDevices = selectedDevices;
                return returnEnv;
        }
	/* 
	catch (cl::Error e) {
                std::cout << e.what() << ": Error code" << e.err() << std::endl;
        }
	*/
        return returnEnv;
}

/*hard code kernel name list*/
struct kernelList {
        //cl::Kernel gpu_try_assign_kernel;
        cl::Kernel gather_kernel;
        cl::Kernel gpu_assign_read_kernel;
        cl::Kernel gpu_count_tempTPM;
        cl::Kernel gpu_count_TPM;
        cl::Kernel gpu_assign_ASE_kernel;
        cl::Kernel gpu_assign_read_ASE_kernel;
        cl::Kernel gpu_count_PSI;
};

/*compile kernels and return*/
std::vector<cl::Kernel> initCompileKernel(std::vector<cl::Device> devices,
                                          cl::Context contexts,
					  std::string programSource
                                          )
{
        std::ifstream programFile(programSource.c_str());
	//std::ifstream programFile("cl_kernel.h");
        std::string programString(std::istreambuf_iterator<char>(programFile),
                                  (std::istreambuf_iterator<char>()));
        cl::Program::Sources source(
            1,
            std::make_pair(programString.c_str(), programString.length() + 1));
        cl::Program program(contexts, source);
        char buildCLFlag[] = "-cl-std=CL1.2 -Werror";
        program.build(devices, buildCLFlag);

        cl::Kernel gather_kernel(program,"gather_kernel");
        cl::Kernel gpu_try_assign_kernel(program, "gpu_try_assign_kernel");
        cl::Kernel gpu_assign_read_kernel(program, "gpu_assign_read_kernel");
        cl::Kernel gpu_count_tempTPM(program, "gpu_count_tempTPM");
        cl::Kernel gpu_count_TPM(program, "gpu_count_TPM");
        cl::Kernel gpu_assign_ASE_kernel(program, "gpu_assign_ASE_kernel");
        cl::Kernel gpu_assign_read_ASE_kernel(program,
                                              "gpu_assign_read_ASE_kernel");
        cl::Kernel gpu_count_PSI(program, "gpu_count_PSI");
        /*build all kernel*/
        std::vector<cl::Kernel> allKernels;
        program.createKernels(&allKernels);

        return allKernels;
}

void checkCLBuild(cl_int &err){
	if (err != CL_SUCCESS){
		printf("Build Fail!\n");
		if (err == CL_INVALID_VALUE){
			printf("CL_INVALID_VALUE\n");
		}
		if (err == CL_INVALID_DEVICE){
			printf("CL_INVALID_DEVICE\n");
		}
		if (err == CL_INVALID_BINARY){
			printf("CL_INVALID_BINARY\n");
		}
		if (err == CL_INVALID_BUILD_OPTIONS){
			printf("CL_INVALID_BUILD_OPTIONS\n");
		}
		if (err == CL_INVALID_OPERATION){
			printf("CL_INVALID_OPERATION\n");
		}
		if (err == CL_COMPILER_NOT_AVAILABLE){
			printf("CL_COMPILER_NOT_AVAILABLE\n");
		}
		if (err == CL_BUILD_PROGRAM_FAILURE){
			printf("CL_BUILD_PROGRAM_FAILURE\n");
		}
		if (err == CL_INVALID_OPERATION){
			printf("CL_INVALID_OPERATION\n");
		}
		if (err == CL_OUT_OF_RESOURCES){
			printf("CL_OUT_OF_RESOURCES\n");
		}
		if (err == CL_OUT_OF_HOST_MEMORY){
			printf("CL_OUT_OF_HOST_MEMORY\n");
		}
		if (err == CL_INVALID_DEVICE){
			printf("CL_INVALID_DEVICE\n");
		}
	}	
}
inline void checkCLKernel(cl_int &err){
	if (err != CL_SUCCESS){
		printf("Kernel Fail!\n");
		if ( err == CL_INVALID_PROGRAM){
			printf("CL_INVALID_PROGRAM\n");
		}
		if ( err == CL_INVALID_PROGRAM_EXECUTABLE){
			printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
		}
		if ( err == CL_INVALID_KERNEL_NAME){
			printf("CL_INVALID_KERNEL_NAME\n");
		}
		if ( err == CL_INVALID_KERNEL_DEFINITION){
			printf("CL_INVALID_KERNEL_DEFINITION\n");
		}
		if ( err == CL_INVALID_VALUE){
			printf("CL_INVALID_VALUE\n");
		}
		if ( err == CL_OUT_OF_RESOURCES){
			printf("CL_OUT_OF_RESOURCES\n");
		}
		if ( err == CL_OUT_OF_HOST_MEMORY){
			printf("CL_OUT_OF_HOST_MEMORY\n");
		}
	}
}
/*compile kernels and return*/
kernelList initCompileKernel_List(std::vector<cl::Device> devices,
                              cl::Context contexts)
{
        //std::ifstream programFile(programSource.c_str());
	std::ifstream programFile("/home/qianjiaqiang/testocl/include/cl_kernel.cl");
        std::string programString(std::istreambuf_iterator<char>(programFile),
                                  (std::istreambuf_iterator<char>()));
        cl::Program::Sources source(
            1,
            std::make_pair(programString.c_str(), programString.length() + 1));
        cl::Program program(contexts, source);
        char buildCLFlag[] = "-cl-std=CL1.2 -Werror";
        cl_int err=program.build(devices, buildCLFlag);
	checkCLBuild(err);
	std::cout << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
        kernelList chipKernel;

        //chipKernel.gpu_try_assign_kernel =
        //   cl::Kernel(program, "gpu_try_assign_kernel");
        chipKernel.gather_kernel=
                cl::Kernel(program,"gather_kernel",&err);
	checkCLKernel(err);
        chipKernel.gpu_assign_read_kernel =
            cl::Kernel(program, "gpu_assign_read_kernel",&err);
	checkCLKernel(err);
        chipKernel.gpu_count_tempTPM = cl::Kernel(program, "gpu_count_tempTPM",&err);
	checkCLKernel(err);
        chipKernel.gpu_count_TPM = cl::Kernel(program, "gpu_count_TPM",&err);
	checkCLKernel(err);
        chipKernel.gpu_assign_ASE_kernel =
            cl::Kernel(program, "gpu_assign_ASE_kernel",&err);
	checkCLKernel(err);
        chipKernel.gpu_assign_read_ASE_kernel =
            cl::Kernel(program, "gpu_assign_read_ASE_kernel",&err);
	checkCLKernel(err);
        chipKernel.gpu_count_PSI = cl::Kernel(program, "gpu_count_PSI",&err);
	checkCLKernel(err);

        return chipKernel;
}

/*
int main()
{
        cl::Context test;
        test = initCppCLSetup();
        return 0;
}
*/
