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
#include <chrono>

#include "xrt.h"
#include "experimental/xrt_kernel.h"

#define NSAMPLES 1024

int main(int argc, char** argv) {
    // Get npoints from argv
    auto NPOINTS = 8;
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
    xrt::bo in_buff[NPOINTS];
    xrt::run run_dm_in[NPOINTS];
    size_t samples_size = sizeof(int16_t) * NSAMPLES * 2; // 32 * 1024

    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // mm2s -> aie
    for (int i = 0; i < NPOINTS; ++ i) {
        // Get reference to the kernels
        dm_in[i] = xrt::kernel(device, uuid, "mm2s:{mm2s_fft_"+ std::to_string(i) +"}");

        // Allocating the input size of sizeIn to MM2S
        in_buff[i] = xrt::bo(device, samples_size, dm_in[i].group_id(0));

        // Write data to compute unit buffers
        in_buff[i].write(&(sample_vector[i * NSAMPLES][0]));

        // Synchronize input buffers data to device global memory
        in_buff[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);

        // Execute the compute units
        run_dm_in[i] = dm_in[i](in_buff[i], nullptr, NSAMPLES/4);
    }

    // aie -> s2mm
    auto dm_out = xrt::kernel(device, uuid, "s2mm:{s2mm_fft_0}");
    auto out_buff = xrt::bo(device, NPOINTS * samples_size, dm_out.group_id(0)); // 32 * 8 * 1024
    auto run_dm_out = dm_out(out_buff, nullptr, NPOINTS * NSAMPLES/4);

    // Wait for kernels to complete
    for (int i = 0; i < NPOINTS; ++ i) {
        run_dm_in[i].wait();
    }
    run_dm_out.wait();

    // Synchronize the output buffer data from the device
    out_buff.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    // Read output buffer data to local buffer
    out_buff.read(fft_result);

    // Stop timer
    auto end_time = std::chrono::high_resolution_clock::now();

    // Output the data
    std::ofstream outfile("DataOutFFT0.txt");
    for (int i = 0; i < NPOINTS * NSAMPLES; i++) {
        outfile << fft_result[i][0] << " " << fft_result[i][1] << std::endl;
    }
    outfile.close();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "TEST PASSED (" << duration.count() << " us)" << std::endl;

    return 0;
}