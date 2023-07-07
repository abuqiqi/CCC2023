#include "fft_kernel.hpp"
#include <cstdio>
#include <aie_api/utils.hpp>
#include <adf.h>

void butterfly(unsigned l, cint16 *x, cint16 *y)
{
    unsigned m = l >> 1, f=N_POINT / l;

    //m是l的一半
    //i是l子序列的第i项
    for (unsigned i = 0; i < m; i += MAX_VEC_LEN)
    {
        vector<cint16, MAX_VEC_LEN> v_omg;
        cint16 *value=(cint16*)omg+f*i;
        for (unsigned j = 0; j < MAX_VEC_LEN; j++,value+=f)
            chess_unroll_loop(MAX_VEC_LEN)
        {
            v_omg[j] = *value;
        }
        for (cint16 *p = x, *p_out = y; p != x + N_POINT; p += l, p_out += l)
        {
            vector<cint16, MAX_VEC_LEN> v_0 = load_v<MAX_VEC_LEN>(p + i);
            vector<cint16, MAX_VEC_LEN> v_1 = load_v<MAX_VEC_LEN>(p + i + m);
            auto acc_t = mul(v_omg, v_1);
            vector<cint16, MAX_VEC_LEN> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
            v_1 = sub(v_0, v_t);
            v_0 = add(v_0, v_t);
            store_v(p_out + i, v_0);
            store_v(p_out + i + m, v_1);
        }
    }
}

void butterfly_256(unsigned l, cint16 *x, cint16 *y)
{
    unsigned m = l >> 1, f=N_POINT / l;

    //m是l的一半
    //i是l子序列的第i项
    for (unsigned i = 0; i < m; i += MAX_VEC_LEN)
    {
        vector<cint16, MAX_VEC_LEN> v_omg;
        cint16 *value=(cint16*)omg_256+i;
        v_omg = load_v<MAX_VEC_LEN>(value);
        // for (unsigned j = 0; j < MAX_VEC_LEN; j++,value+=f)
        //     chess_unroll_loop(MAX_VEC_LEN)
        // {
        //     v_omg[j] = *value;
        // }
        for (cint16 *p = x, *p_out = y; p != x + N_POINT; p += l, p_out += l)
        {
            vector<cint16, MAX_VEC_LEN> v_0 = load_v<MAX_VEC_LEN>(p + i);
            vector<cint16, MAX_VEC_LEN> v_1 = load_v<MAX_VEC_LEN>(p + i + m);
            auto acc_t = mul(v_omg, v_1);
            vector<cint16, MAX_VEC_LEN> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
            v_1 = sub(v_0, v_t);
            v_0 = add(v_0, v_t);
            store_v(p_out + i, v_0);
            store_v(p_out + i + m, v_1);
        }
    }
}

void butterfly_512(unsigned l, cint16 *x, cint16 *y)
{
    unsigned m = l >> 1, f=N_POINT / l;

    //m是l的一半
    //i是l子序列的第i项
    for (unsigned i = 0; i < m; i += MAX_VEC_LEN)
    {
        vector<cint16, MAX_VEC_LEN> v_omg;
        vector<cint16, 16> v_tmp;
        cint16 *value=(cint16*)omg+f*i;
        v_tmp = load_v<16>(value);
        v_omg = filter_even(v_tmp);
        // for (unsigned j = 0; j < MAX_VEC_LEN; j++,value+=f)
        //     chess_unroll_loop(MAX_VEC_LEN)
        // {
        //     v_omg[j] = *value;
        // }
        for (cint16 *p = x, *p_out = y; p != x + N_POINT; p += l, p_out += l)
        {
            vector<cint16, MAX_VEC_LEN> v_0 = load_v<MAX_VEC_LEN>(p + i);
            vector<cint16, MAX_VEC_LEN> v_1 = load_v<MAX_VEC_LEN>(p + i + m);
            auto acc_t = mul(v_omg, v_1);
            vector<cint16, MAX_VEC_LEN> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
            v_1 = sub(v_0, v_t);
            v_0 = add(v_0, v_t);
            store_v(p_out + i, v_0);
            store_v(p_out + i + m, v_1);
        }
    }
}

void butterfly_1024(unsigned l, cint16 *x, cint16 *y)
{
    unsigned m = l >> 1, f=N_POINT / l;

    //m是l的一半
    //i是l子序列的第i项
    for (unsigned i = 0; i < m; i += MAX_VEC_LEN)
    {
        vector<cint16, MAX_VEC_LEN> v_omg;
        cint16 *value=(cint16*)omg+i;
        v_omg = load_v<MAX_VEC_LEN>(value);
        // for (unsigned j = 0; j < MAX_VEC_LEN; j++,value+=f)
        //     chess_unroll_loop(MAX_VEC_LEN)
        // {
        //     v_omg[j] = *value;
        // }
        vector<cint16, MAX_VEC_LEN> v_0 = load_v<MAX_VEC_LEN>(x + i);
        vector<cint16, MAX_VEC_LEN> v_1 = load_v<MAX_VEC_LEN>(x + i + m);
        auto acc_t = mul(v_omg, v_1);
        vector<cint16, MAX_VEC_LEN> v_t = acc_t.to_vector<cint16>(OMG_SHIFT);
        v_1 = sub(v_0, v_t);
        v_0 = add(v_0, v_t);
        store_v(y + i, v_0);
        store_v(y + i + m, v_1);
    }
}

// void mul_twiddle_factor(cint16 *x,cint16 *y,cint16 *tf)    
// {
//     auto iter_tf=begin_vector<MAX_VEC_LEN>(tf);
//     auto iter_x=begin_vector<MAX_VEC_LEN>(x);
//     auto iter_y=begin_vector<MAX_VEC_LEN>(y);
//     // for (cint16 *p_out=y;p_out!=y+N_POINT;p_out+=MAX_VEC_LEN)
//     for (unsigned i=0;i<N_POINT;i+=MAX_VEC_LEN)
//     {
//         auto res=mul(*iter_x++,*iter_tf++);
//         // store_v(p_out,res.to_vector<cint16>(TF_SHIFT));
//         // window_writeincr(y,res.to_vector<cint16>(TF_SHIFT));
//         *iter_y++=res.to_vector<cint16>(TF_SHIFT);
//     }
// }

void mul_twiddle_factor(cint16 *x,cint16 *y,cint16 *tf)    
{
    auto iter_tf=begin_vector<32>(tf);
    auto iter_x=begin_vector<32>(x);
    auto iter_y=begin_vector<32>(y);
    // for (cint16 *p_out=y;p_out!=y+N_POINT;p_out+=MAX_VEC_LEN)
    for (unsigned i=0;i<N_POINT;i+=32)
    {
        auto res=mul(*iter_x++,*iter_tf++);
        // store_v(p_out,res.to_vector<cint16>(TF_SHIFT));
        // window_writeincr(y,res.to_vector<cint16>(TF_SHIFT));
        *iter_y++=res.to_vector<cint16>(TF_SHIFT);
    }
}

template<unsigned id>
void radix2_dit(input_window<cint16> *x_in, output_window<cint16> *y_out)
{
    cint16 *x = (cint16 *)x_in->ptr;
    cint16 *y = (cint16 *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    aie::tile tile = aie::tile::current();
    printf("before shuffle: %llu\n", tile.cycles());

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
        chess_unroll_loop(8)
    {
        int16 a=swap2[i],b=swap2[i+1];
        // cint16 tmp = x[a];
        // x[a] = x[b];
        // x[b] = tmp;
        y[a] = x[b];
        y[b] = x[a];
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
        chess_unroll_loop(8)
    {
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        // cint16 tmp = x[i3];
        // x[i3] = x[i0];
        // x[i0] = x[i1];
        // x[i1] = x[i2];
        // x[i2] = tmp;
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

    printf("shuffled: %llu\n", tile.cycles());

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

    // auto iterin=begin_vector<MAX_VEC_LEN>(x);
    // auto iterout=begin_vector<MAX_VEC_LEN>(y);
    // for (cint16* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN)
    // {
    //     auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
    //     auto 
    //     auto m=mul(*iter++,*p);

    //     for (unsigned i=1;i<MAX_VEC_LEN;i++){
    //         m=mac(m,*iter++,*(p+i));
    //     }
    //     // store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
    //     *iterout++=m.to_vector<cint16>(MAT_OMG_SHIFT);
    // }

    // auto iterx=begin_vector<4>(x);
    // for (cint16* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN ){
    //     vector<cint16,32> v_omg=load_v<32>(omg8_0);
    //     auto m=sliding_mul_ops<8,4,1,1,4,cint16,cint16,cacc48>::mul(*iterx++,0,v_omg,0);
    //     v_omg=load_v<32>(omg8_1);
    //     m=sliding_mul_ops<8,4,1,1,4,cint16,cint16,cacc48>::mac(m,*iterx++,0,v_omg,0);
    //     store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
    // }

    printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    

    if (id==0){
        // for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT >> 3; l <<= 1)
        // {
        //     butterfly(l, y, y);
        //     printf("btf l=%u: %llu\n", l, tile.cycles());
        // }
        butterfly(16, x, y);
        printf("id:0 btf l=16: %llu\n", tile.cycles());
        butterfly(32, y, x);
        printf("id:0 btf l=32: %llu\n", tile.cycles());
        butterfly(64, x, y);
        printf("id:0 btf l=64: %llu\n", tile.cycles());
        butterfly(128, y, x);
        printf("id:0 btf l=128: %llu\n", tile.cycles());
        butterfly_256(256, x, y);
        printf("id:0 btf l=256: %llu\n", tile.cycles());
        butterfly_512(512, y, x);
        printf("id:0 btf l=512: %llu\n", tile.cycles());
        butterfly_1024(N_POINT, x, y);
        printf("id:0 btf l=%u: %llu\n", N_POINT, tile.cycles());
    }
    else {
        // for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT >> 3; l <<= 1)
        // {
        //     butterfly(l, y, y);
        //     printf("btf l=%u: %llu\n", l, tile.cycles());
        // }
        butterfly(16, x, y);
        printf("btf l=16: %llu\n", tile.cycles());
        butterfly(32, y, x);
        printf("btf l=32: %llu\n", tile.cycles());
        butterfly(64, x, y);
        printf("btf l=64: %llu\n", tile.cycles());
        butterfly(128, y, x);
        printf("btf l=128: %llu\n", tile.cycles());
        butterfly_256(256, x, y);
        printf("btf l=256: %llu\n", tile.cycles());
        butterfly_512(512, y, x);
        printf("btf l=512: %llu\n", tile.cycles());
        butterfly_1024(N_POINT, x, y);
        printf("btf l=%u: %llu\n", N_POINT, tile.cycles());
        switch (id)
        {
        case 1:
            mul_twiddle_factor(y,y,(cint16*)tf1);
            break;
        case 2:
            mul_twiddle_factor(y,y,(cint16*)tf2);
            break;
        case 3:
            mul_twiddle_factor(y,y,(cint16*)tf3);
            break;
        }
    }

    printf("dit: %llu\n", tile.cycles());

    return;
}