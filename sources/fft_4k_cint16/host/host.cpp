// Copyright (C) 2023 Advanced Micro Devices, Inc
//
// SPDX-License-Identifier: MIT

#include <assert.h>
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstdint>

#include "xrt.h"
#include "experimental/xrt_kernel.h"

#define NSAMPLES 1024

int main(int argc, char** argv) {
    // Get npoints from argv
    auto NPOINTS = 4;
    if ( argc == 2 ) {
        NPOINTS = std::stoi(argv[1]);
    }
    std::cout << "Load the point size " << NPOINTS << "*" << NSAMPLES << std::endl;

    // Get device index and download xclbin
    std::cout << "Open the device" << std::endl;
    auto device = xrt::device(0);
    std::string binaryFile = "./fft.xclbin";
    std::cout << "Load the xclbin " << binaryFile << std::endl;
    auto uuid = device.load_xclbin(binaryFile);

    // Read generated data
    auto *sample_vector = new int16_t [NPOINTS * NSAMPLES][2];
    auto *fft_result = new int16_t [NPOINTS * NSAMPLES][2];

    std::ifstream infile("DataInFFT0.txt");
    for (int i = 0; i < NPOINTS * NSAMPLES; i++) {
        infile >> sample_vector[i][0] >> sample_vector[i][1];
    }
    infile.close();

    // Create datamover objects, in/out buffer
    xrt::kernel dm_in[NPOINTS];
    xrt::kernel dm_out[NPOINTS];
    xrt::bo in_buff[NPOINTS];
    xrt::bo out_buff[NPOINTS];
    xrt::run run_dm_in[NPOINTS];
    xrt::run run_dm_out[NPOINTS];
    size_t samples_size = sizeof(int16_t) * NSAMPLES * 2;

    for (int i = 0; i < NPOINTS; ++ i) {
        // Get reference to the kernels
        std::cout << "Get references to datamovers compute units" << std::endl;
        dm_in[i] = xrt::kernel(device, uuid, "mm2s:{mm2s_fft_"+ std::to_string(i) +"}");
        dm_out[i] = xrt::kernel(device, uuid, "s2mm:{s2mm_fft_"+ std::to_string(i) +"}");

        // Allocating the input size of sizeIn to MM2S
        std::cout << "Allocate Buffer in Global Memory" << std::endl;
        in_buff[i] = xrt::bo(device, samples_size, dm_in[i].group_id(0));
        out_buff[i] = xrt::bo(device, samples_size, dm_out[i].group_id(0));

        // Write data to compute unit buffers
        in_buff[i].write(&(sample_vector[i * NSAMPLES][0]));

        // Synchronize input buffers data to device global memory
        in_buff[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);

        // Execute the compute units
        run_dm_in[i] = dm_in[i](in_buff[i], nullptr, NSAMPLES/4); // mm2s -> aie
        run_dm_out[i] = dm_out[i](out_buff[i], nullptr, NSAMPLES/4); // aie -> s2mm
        std::cout << "Kernels started" << std::endl;
    }

    // Wait for kernels to complete
    for (int i = 0; i < NPOINTS; ++ i) {
        run_dm_in[i].wait();
        run_dm_out[i].wait();
    }

    for (int i = 0; i < NPOINTS; ++ i) {
        // Synchronize the output buffer data from the device
        out_buff[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);

        // Read output buffer data to local buffer
        out_buff[i].read(&(fft_result[i * NSAMPLES][0]));
    }

    // Output the data
    std::ofstream outfile("DataOutFFT0.txt");
    for (int i = 0; i < NPOINTS * NSAMPLES; i++) {
        outfile << fft_result[i][0] << " " << fft_result[i][1] << std::endl;
    }
    outfile.close();

    std::cout << "TEST PASSED" << std::endl;

    return 0;

    // auto mm2s_0 = xrt::kernel(device, uuid, "mm2s:{mm2s_fft_0}");
    // auto s2mm_0 = xrt::kernel(device, uuid, "s2mm:{s2mm_fft_0}");

    // auto in_buff_1 = xrt::bo(device, samples_size, mm2s_0.group_id(0));
    // auto out_buff_1 = xrt::bo(device, samples_size, s2mm_0.group_id(0));

    // in_buff_1.write(sample_vector); // infile -> sample_vector -> in_buff_1
    
    // in_buff_1.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // auto run_mm2s_0 = mm2s_0(in_buff_1, nullptr, NSAMPLES/4);
    // auto run_s2mm_0 = s2mm_0(out_buff_1, nullptr, NSAMPLES/4);

    // run_mm2s_0.wait();
    // run_s2mm_0.wait();

}