#include "fft_kernel.hpp"
#include <aie_api/aie.hpp>
#include <cstdio>
#include <aie_api/utils.hpp>

#define PI_ELEMENTARY 3.14159265359

using namespace aie;

cfloat omega(int k){
    return cfloat(cos((float)2.0 * (float)PI_ELEMENTARY * k / N_POINT), (-1)*sin((float)2.0 * (float)PI_ELEMENTARY * k / N_POINT));
}

void radix2_dit(input_stream<cfloat> * x_in,output_stream<cfloat> * y_out){
    cfloat x[N_POINT];
    for (int i=0;i<N_POINT;i++) x[i]=readincr(x_in);

    int lim = 0;
    while((1 << lim) < N_POINT) lim++;
    printf("lim: %d\n",lim);
    for(int i = 0; i < N_POINT; i++){
    	int t = 0;
    	for(int j = 0; j < lim; j++)
    	    if((i >> j) & 1) t |= (1 << (lim - j - 1));
    	if(i < t){
            cfloat tmp=x[i];
            x[i]=x[t];
            x[t]=tmp;
        }
    }
    // printf("swap\n");
    // for (int i=0;i<N_POINT;i++) printf("%f %f ",x[i].real,x[i].imag);
    // printf("\n");
    for(int l = 2; l <= N_POINT; l *= 2){
    	int m = l / 2;
        for (cfloat *p=x;p!=x+N_POINT;p+=l)
            for(int i = 0; i < m; i++){
                cfloat t = (cfloat)omega(N_POINT / l * i) * p[i + m];
                // printf("t=%f %f\n",t.real,t.imag);
                // cfloat omg = omega(N_POINT / l * i);
                // printf("omg=%f %f\n",omg.real,omg.imag);
                p[i + m] = p[i] - t;
                p[i] += t;
            }
        // printf("l = %d\n",l);
        // for (int i=0;i<N_POINT;i++) printf("%f %f ",x[i].real,x[i].imag);
        // printf("\n");
    }

    for (int i=0;i<N_POINT;i++) writeincr(y_out,x[i]);
    // for (int i=0;i<N_POINT;i+=4){
    //     vector<cfloat,4> v=load_v<4>((cfloat*)(x+i));
    //     writeincr(y_out,v);
    // }
    return;
}
