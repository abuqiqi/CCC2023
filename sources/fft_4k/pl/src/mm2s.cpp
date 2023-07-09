// Copyright (C) 2023 Advanced Micro Devices, Inc
//
// SPDX-License-Identifier: MIT

#include "config.hpp"

extern "C" {

void mm2s(ap_int<DWIDTH>* mem, hls::stream<data >& s, int size) {
data_mover:
    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1 // pipeline

        // #pragma HLS loop_tripcount min=512 max=524288 // cint16 - 32x64/4 to 1024x2048/4

        data x;
        x.data = mem[i];
        x.keep_all();
        s.write(x);
    }
}
}