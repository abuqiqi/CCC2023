#pragma once

#include <adf.h>
#include "fft.hpp"

using namespace adf;

class main_graph: public graph{
private:
    aie_fft_graph fft;
public:
    input_plio in;
    output_plio out;
    
    main_graph(){
        in=input_plio::create("DataInFFT0",plio_128_bits,"data/DataInFFT0.txt");
        out=output_plio::create("DataOutFFT0",plio_128_bits,"data/DataOutFFT0.txt");
        connect<>(in.out[0],fft.in);
        connect<>(fft.out,out.in[0]);
    }
};

