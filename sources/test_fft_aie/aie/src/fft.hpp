#pragma once

#include <adf.h>
#include <fft_kernel.hpp>

using namespace adf;

class aie_fft_graph : public graph {
private:
    kernel fft_kernel;
public:
    port<input> in;
    port<output> out;

    aie_fft_graph(){
        fft_kernel=kernel::create(radix2_dit);

        connect<>(in,fft_kernel.in[0]);
        // connect<>(tw,fft_kernel.in[1]);
        connect<>(fft_kernel.out[0],out);

        source(fft_kernel)="fft_kernel.cpp";

        runtime<ratio>(fft_kernel)=1;
    }
};