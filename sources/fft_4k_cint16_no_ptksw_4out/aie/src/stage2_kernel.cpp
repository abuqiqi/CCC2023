#include "stage2_kernel.hpp"
#include <aie_api/utils.hpp>
#include <cstdio>

void fft_stage2(input_window<cint16> *x_in0,input_window<cint16> *x_in1,input_window<cint16> *x_in2,input_window<cint16> *x_in3,
                output_window<cint16> *y_out0,output_window<cint16> *y_out1,output_window<cint16> *y_out2,output_window<cint16> *y_out3)
{
    aie::tile tile = aie::tile::current();
    printf("before stage2: %llu\n", tile.cycles());

    cint16 *x0=(cint16*)x_in0->ptr;
    cint16 *x1=(cint16*)x_in1->ptr;
    cint16 *x2=(cint16*)x_in2->ptr;
    cint16 *x3=(cint16*)x_in3->ptr;
    cint16 *y0=(cint16*)y_out0->ptr;
    cint16 *y1=(cint16*)y_out1->ptr;
    cint16 *y2=(cint16*)y_out2->ptr;
    cint16 *y3=(cint16*)y_out3->ptr;
    for (unsigned i=0;i<N_POINT;i++){
        auto iter=begin_vector<4>(omg4);
        auto m=mul(*iter++,x0[i]);
        m=mac(m,*iter++,x1[i]);
        m=mac(m,*iter++,x2[i]);
        m=mac(m,*iter++,x3[i]);
        vector<cint16,4> v=m.to_vector<cint16>(0);
        *y0++=v[0];
        *y1++=v[1];
        *y2++=v[2];
        *y3++=v[3];
    }

    printf("stage2: %llu\n", tile.cycles());
}