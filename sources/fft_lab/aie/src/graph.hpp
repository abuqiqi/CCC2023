// Copyright (C) 2023 Advanced Micro Devices, Inc
//
// SPDX-License-Identifier: MIT

#pragma once

#include "adf.h"
#include "fft.hpp"

class SpectrumGraph: public graph {
    private:
        FFT1d_graph fft;

    public:
        adf::port<input> in;
        adf::port<output> out;
        SpectrumGraph() {
            connect<stream>(in, fft.in);
            // connect< window<FFT_POINT_SIZE * 4> > (fir.out,  fft.in);
            connect<stream>(fft.out, out);
        }
};


class DSPLibGraph: public graph {
    private:
        FFT1d_graph fft;
        // SpectrumGraph spectrum;
        
    public:
        adf::input_plio s_in[1];
        adf::output_plio s_out[1];
        DSPLibGraph() {
            // FFT connections
            s_in[0] = input_plio::create("DataInFFT0", adf::plio_128_bits, "data/DataInFFT0.txt");
            s_out[0] = output_plio::create("DataOutFFT0", adf::plio_128_bits, "DataOutFFT0.txt");
            adf::connect<stream> (s_in[0].out[0],  fft.in);
            adf::connect<stream> (fft.out, s_out[0].in[0]);
        }
};
