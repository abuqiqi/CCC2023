#include "fft_kernel.hpp"
#include <cstdio>
#include <aie_api/utils.hpp>
#include <adf.h>

// #define PI_ELEMENTARY 3.14159265359

// cfloat omega(unsigned k){
//     return cfloat(cos((float)2.0 * (float)PI_ELEMENTARY * k / N_POINT), (-1)*sin((float)2.0 * (float)PI_ELEMENTARY * k / N_POINT));
// }

void butterfly(unsigned l, cfloat *x, cfloat *y)
{
    // aie::tile tile = aie::tile::current();
    unsigned m = l >> 1, f=N_POINT / l;
    for (unsigned i = 0; i < m; i += MAX_VEC_LEN)
        // chess_prepare_for_pipelining
        // chess_loop_range(1,)
    {
        vector<cfloat, MAX_VEC_LEN> v_omg;
        cfloat *value=(cfloat*)omg+f*i;
        for (unsigned j = 0; j < MAX_VEC_LEN; j++,value+=f)
            chess_unroll_loop(MAX_VEC_LEN)
        {
            v_omg[j] = *value;
        }
        // printf("get omg: %llu\n", tile.cycles());
        for (cfloat *p = x, *p_out = y; p != x + N_POINT; p += l, p_out += l)
            // chess_prepare_for_pipelining
            // chess_loop_range(1,)
        {
            vector<cfloat, MAX_VEC_LEN> v_0 = load_v<MAX_VEC_LEN>(p + i);
            vector<cfloat, MAX_VEC_LEN> v_1 = load_v<MAX_VEC_LEN>(p + i + m);
            auto acc_t = mul(v_omg, v_1);
            vector<cfloat, MAX_VEC_LEN> v_t = acc_t.to_vector<cfloat>(0);
            v_1 = sub(v_0, v_t);
            v_0 = add(v_0, v_t);
            store_v(p_out + i, v_0);
            store_v(p_out + i + m, v_1);
        }
        // printf("compute: %llu\n", tile.cycles());
    }
}

void radix2_dit(input_window<cfloat> *x_in, output_window<cfloat> *y_out)
{
    cfloat *x = (cfloat *)x_in->ptr;
    cfloat *y = (cfloat *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    // for(unsigned i = 0; i < N_POINT; i++){
    // 	unsigned t = 0;
    // 	for(unsigned j = 0; j < lim; j++)
    // 	    if((i >> j) & 1) t |= (1 << (lim - j - 1));
    // 	if(i < t){
    //         cfloat tmp=x[i];
    //         x[i]=x[t];
    //         x[t]=tmp;
    //     }
    // }

    aie::tile tile = aie::tile::current();
    printf("before shuffle: %llu\n", tile.cycles());

    // for (unsigned i = 0; i < N_SWAP; i += 2)
    //     chess_prepare_for_pipelining   
    //     chess_loop_count(N_SWAP/2)
    // {
    //     cfloat tmp = x[swap[i]];
    //     x[swap[i]] = x[swap[i + 1]];
    //     x[swap[i + 1]] = tmp;
    // }
    // for (cfloat * p=x;p!=x+N_POINT;p+=4)
    //     chess_prepare_for_pipelining  
    //     chess_loop_count(N_SWAP/4)
    // {
    //     cfloat tmp = *(p+1);
    //     *(p+1) = *(p+2);
    //     *(p+2)=tmp;
    // }
    // for (cfloat *p=x;p!=x+N_POINT;p+=8)
    //     chess_prepare_for_pipelining   
    //     chess_loop_count(N_SWAP/8)
    // {
    //     cfloat tmp=*(p+4);
    //     *(p+4)=*(p+2);
    //     *(p+2)=*(p+1);
    //     *(p+1)=tmp;
    //     tmp=*(p+6);
    //     *(p+6)=*(p+3);
    //     *(p+3)=*(p+5);
    //     *(p+5)=tmp;
    // }

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
        // chess_prepare_for_pipelining   
        // chess_loop_count(N_SWAP_2/2)
        // chess_unroll_loop(2)
    {
        // cfloat tmp = x[swap2[i]];
        // x[swap2[i]] = x[swap2[i + 1]];
        // x[swap2[i + 1]] = tmp;
        int16 a=swap2[i],b=swap2[i+1];
        cfloat tmp = x[a];
        x[a] = x[b];
        x[b] = tmp;
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
        // chess_prepare_for_pipelining   
        // chess_loop_count(N_SWAP_4/4)
        // chess_unroll_loop(8)
    {
        // cfloat tmp = x[swap4[i+3]];
        // x[swap4[i+3]] = x[swap4[i]];
        // x[swap4[i]] = x[swap4[i+1]];
        // x[swap4[i+1]] = x[swap4[i+2]];
        // x[swap4[i + 2]] = tmp;
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        cfloat tmp = x[i3];
        x[i3] = x[i0];
        x[i0] = x[i1];
        x[i1] = x[i2];
        x[i2] = tmp;
    }

    // for (unsigned i=0;i<N_POINT;i++) printf("(%f,%f),",x[i].real,x[i].imag);
    // printf("\n");

    printf("shuffled: %llu\n", tile.cycles());

    // ----------------------------------dit----------------------------------

    // for(unsigned l = 2; l <= N_POINT; l *= 2){
    // 	unsigned m = l / 2;
    //     for (cfloat *p=x;p!=x+N_POINT;p+=l)
    //         for(unsigned i = 0; i < m; i++){
    //             // cfloat t = omg[N_POINT / l * i] * p[i + m];
    //             cfloat t = omega(N_POINT / l * i) * p[i + m];
    //             // printf("t=%f %f\n",t.real,t.imag);
    //             // cfloat omg = omega(N_POINT / l * i);
    //             // printf("omg=%f %f\n",omg.real,omg.imag);
    //             p[i + m] = p[i] - t;
    //             p[i] += t;
    //         }
    //     // printf("l = %d\n",l);
    //     // for (unsigned i=0;i<N_POINT;i++) printf("%f %f ",x[i].real,x[i].imag);
    //     // printf("\n");
    // }

    for (cfloat* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN)
        // chess_prepare_for_pipelining   
        // chess_loop_count(N_POINT/MAX_VEC_LEN)
        // chess_unroll_loop(2)
    {
        auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
        auto m=mul(*iter++,*p);

        for (unsigned i=1;i<MAX_VEC_LEN;i++){
            m=mac(m,*iter++,*(p+i));
        }
        store_v(p,m.to_vector<cfloat>(0));
    }
    printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT >> 1; l <<= 1)
        // chess_prepare_for_pipelining   
        // chess_loop_count(7)
    {
        butterfly(l, x, x);
        printf("btf l=%u: %llu\n", l, tile.cycles());
    }

    // !!! last round
    butterfly(N_POINT, x, y);

    printf("dit: %llu\n", tile.cycles());

    return;
}