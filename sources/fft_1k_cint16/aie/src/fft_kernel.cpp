#include "fft_kernel.hpp"
#include <cstdio>
#include <aie_api/utils.hpp>
#include <adf.h>

void butterfly(unsigned l, cint16 *x, cint16 *y, cint16 *omg)
{
    unsigned m = l >> 1;
	auto iteromg=begin_vector<32>(omg);
    for (unsigned i = 0; i < m; i += 32)
    {
        for (cint16 *p = x, *p_out = y; p != x + N_POINT; p += l, p_out += l)
        {
            vector<cint16, 32> v_0 = load_v<32>(p + i);
            vector<cint16, 32> v_1 = load_v<32>(p + i + m);
            // if (l==1024 && i==128) printf("before: x[139] = %hd+%hdi, x[139+512] = %hd+%hdi\n",((cint16)v_0[11]).real,((cint16)v_0[11]).imag,((cint16)v_1[11]).real,((cint16)v_1[11]).imag);
            // if (l==1024 && i==128) printf("omg = %hd+%hdi\n",((cint16)((*iteromg)[11])).real,((cint16)((*iteromg)[11])).imag);
            auto acc_t = mul(*iteromg, v_1);
            // if (l==1024 && i==128) print(acc_t,true,"acc_t: ");
            vector<cint16, 32> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
            // if (l==1024 && i==128) print(v_t,true,"v_t: ");
            // if (l==1024 && i==128) printf("x[139+512]*omg = %x+%hdi\n",((cint16)v_t[11]).real,((cint16)v_t[11]).imag);
            v_1 = sub(v_0, v_t);
            v_0 = add(v_0, v_t);
            // if (l==1024 && i==128) printf("after: x[139] = %hd+%hdi, x[139+512] = %hd+%hdi\n",((cint16)v_0[11]).real,((cint16)v_0[11]).imag,((cint16)v_1[11]).real,((cint16)v_1[11]).imag);
            store_v(p_out + i, v_0);
            store_v(p_out + i + m, v_1);

            // accum<cacc48,32> acc_a;
            // acc_a.from_vector(v_0,0);
            // store_v(p_out + i + m, sub(acc_a,v_t).to_vector<cint16>(0));
            // store_v(p_out + i, add(acc_a,v_t).to_vector<cint16>(0));

            // accum<cacc48,32> acc_a;
            // acc_a.from_vector(v_0,OMG_SHIFT);
            // store_v(p_out + i + m, msc(acc_a,*iteromg, v_1).to_vector<cint16>(OMG_SHIFT));
            // store_v(p_out + i, mac(acc_a,*iteromg, v_1).to_vector<cint16>(OMG_SHIFT));
        }
	    iteromg++;
    }
}

void butterfly_16(cint16 *x, cint16 *y)
{
	vector<cint16, 8> v_omg = load_v<8>(omg_16);
    for (cint16 *p = x, *p_out = y; p != x + N_POINT; p += 16, p_out += 16)
    {
        vector<cint16, 8> v_0 = load_v<8>(p);
        vector<cint16, 8> v_1 = load_v<8>(p + 8);
        auto acc_t = mul(v_omg, v_1);
        vector<cint16, 8> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
        v_1 = sub(v_0, v_t);
        v_0 = add(v_0, v_t); 
        store_v(p_out, v_0);
        store_v(p_out + 8, v_1);
    }
}

void butterfly_32(cint16 *x, cint16 *y)
{
	vector<cint16, 16> v_omg = load_v<16>(omg_32);
    for (cint16 *p = x, *p_out = y; p != x + N_POINT; p += 32, p_out += 32)
    {
        vector<cint16, 16> v_0 = load_v<16>(p);
        vector<cint16, 16> v_1 = load_v<16>(p + 16);
        auto acc_t = mul(v_omg, v_1);
        vector<cint16, 16> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
        v_1 = sub(v_0, v_t);
        v_0 = add(v_0, v_t);
        store_v(p_out, v_0);
        store_v(p_out + 16, v_1);
    }
}

void radix2_dit(input_window<cint16> *x_in, output_window<cint16> *y_out)
{
    cint16 *x = (cint16 *)x_in->ptr;
    cint16 *y = (cint16 *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    // aie::tile tile = aie::tile::current();
    // // printf("before shuffle: %llu\n", tile.cycles());

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
        chess_unroll_loop(8)
    {
        int16 a=swap2[i],b=swap2[i+1];
        y[a] = x[b];
        y[b] = x[a];
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
        chess_unroll_loop(8)
    {
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        y[i3] = x[i0];
        y[i0] = x[i1];
        y[i1] = x[i2];
        y[i2] = x[i3];
    }
    //没换的要拷贝过去
    y[0] = x[0];
    y[48] = x[48];
    y[72] = x[72];
    y[120] = x[120];
    y[258] = x[258];
    y[306] = x[306];
    y[330] = x[330];
    y[378] = x[378];
    y[645] = x[645];
    y[693] = x[693];
    y[717] = x[717];
    y[765] = x[765];
    y[903] = x[903];
    y[951] = x[951];
    y[975] = x[975];
    y[1023] = x[1023];

    // printf("shuffled: %llu\n", tile.cycles());

    // ----------------------------------dit----------------------------------

    auto iterout=begin_vector<MAX_VEC_LEN>(x);
    for (cint16* p=y;p!=y+N_POINT;p+=MAX_VEC_LEN)
    {
        auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
        auto m=mul(*iter++,*p);

        for (unsigned i=1;i<MAX_VEC_LEN;i++){
            m=mac(m,*iter++,*(p+i));
        }
        // store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
        *iterout++=m.to_vector<cint16>(MAT_OMG_SHIFT);
    }

    // printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    
    butterfly_16(x, y);
    // printf("btf l=16: %llu\n", tile.cycles());
    butterfly_32(y, x);
    // printf("btf l=32: %llu\n", tile.cycles());
    butterfly(64, x, y,omg_64);
    // printf("btf l=64: %llu\n", tile.cycles());
    butterfly(128, y, x,omg_128);
    // printf("btf l=128: %llu\n", tile.cycles());
    butterfly(256, x, y,omg_256);
    // printf("btf l=256: %llu\n", tile.cycles());
    butterfly(512, y, x,omg_512);
    // printf("btf l=512: %llu\n", tile.cycles());
    butterfly(1024, x, y,omg_1024);
    // printf("btf l=1024: %llu\n", tile.cycles());

    // printf("dit: %llu\n", tile.cycles());

    return;
}