#include "fft_kernel.hpp"
#include <cstdio>
#include <aie_api/utils.hpp>
#include <adf.h>

void butterfly(unsigned l, cint16 *x, cint16 *y)
{
    unsigned m = l >> 1, f=N_POINT / l;
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

void mul_twiddle_factor(cint16 *x,output_window<cint16> *y,cint16 *tf)    
{
    auto iter_tf=begin_vector<MAX_VEC_LEN>(tf);
    auto iter_x=begin_vector<MAX_VEC_LEN>(x);
    // for (cint16 *p_out=y;p_out!=y+N_POINT;p_out+=MAX_VEC_LEN)
    for (unsigned i=0;i<N_POINT;i+=MAX_VEC_LEN)
    {
        auto res=mul(*iter_x++,*iter_tf++);
        // store_v(p_out,res.to_vector<cint16>(TF_SHIFT));
        window_writeincr(y,res.to_vector<cint16>(TF_SHIFT));
    }
}

void radix2_dit_0(input_window<cint16> *x_in, output_window<cint16> *y_out)
{
    cint16 *x = (cint16 *)x_in->ptr;
    cint16 *y = (cint16 *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    aie::tile tile = aie::tile::current();
    printf("before shuffle: %llu\n", tile.cycles());

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
    {
        int16 a=swap2[i],b=swap2[i+1];
        cint16 tmp = x[a];
        x[a] = x[b];
        x[b] = tmp;
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
    {
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        cint16 tmp = x[i3];
        x[i3] = x[i0];
        x[i0] = x[i1];
        x[i1] = x[i2];
        x[i2] = tmp;
    }

    printf("shuffled: %llu\n", tile.cycles());

    // ----------------------------------dit----------------------------------

    for (cint16* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN)
    {
        auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
        auto m=mul(*iter++,*p);

        for (unsigned i=1;i<MAX_VEC_LEN;i++){
            m=mac(m,*iter++,*(p+i));
        }
        store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
    }
    printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT >> 1; l <<= 1)
    {
        butterfly(l, x, x);
        printf("btf l=%u: %llu\n", l, tile.cycles());
    }

    butterfly(N_POINT, x, y);

    printf("dit: %llu\n", tile.cycles());

    return;
}

void radix2_dit_1(input_window<cint16> *x_in, output_window<cint16> *y_out)
{
    cint16 *x = (cint16 *)x_in->ptr;
    cint16 *y = (cint16 *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    aie::tile tile = aie::tile::current();
    printf("before shuffle: %llu\n", tile.cycles());

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
    {
        int16 a=swap2[i],b=swap2[i+1];
        cint16 tmp = x[a];
        x[a] = x[b];
        x[b] = tmp;
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
    {
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        cint16 tmp = x[i3];
        x[i3] = x[i0];
        x[i0] = x[i1];
        x[i1] = x[i2];
        x[i2] = tmp;
    }

    printf("shuffled: %llu\n", tile.cycles());

    // ----------------------------------dit----------------------------------

    for (cint16* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN)
    {
        auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
        auto m=mul(*iter++,*p);

        for (unsigned i=1;i<MAX_VEC_LEN;i++){
            m=mac(m,*iter++,*(p+i));
        }
        store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
    }
    printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT ; l <<= 1)
    {
        butterfly(l, x, x);
        printf("btf l=%u: %llu\n", l, tile.cycles());
    }

    mul_twiddle_factor(x,y_out,(cint16*)tf1);

    printf("dit: %llu\n", tile.cycles());

    return;
}

void radix2_dit_2(input_window<cint16> *x_in, output_window<cint16> *y_out)
{
    cint16 *x = (cint16 *)x_in->ptr;
    cint16 *y = (cint16 *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    aie::tile tile = aie::tile::current();
    printf("before shuffle: %llu\n", tile.cycles());

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
    {
        int16 a=swap2[i],b=swap2[i+1];
        cint16 tmp = x[a];
        x[a] = x[b];
        x[b] = tmp;
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
    {
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        cint16 tmp = x[i3];
        x[i3] = x[i0];
        x[i0] = x[i1];
        x[i1] = x[i2];
        x[i2] = tmp;
    }

    printf("shuffled: %llu\n", tile.cycles());

    // ----------------------------------dit----------------------------------

    for (cint16* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN)
    {
        auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
        auto m=mul(*iter++,*p);

        for (unsigned i=1;i<MAX_VEC_LEN;i++){
            m=mac(m,*iter++,*(p+i));
        }
        store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
    }
    printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT; l <<= 1)
    {
        butterfly(l, x, x);
        printf("btf l=%u: %llu\n", l, tile.cycles());
    }

    mul_twiddle_factor(x,y_out,(cint16*)tf2);

    printf("dit: %llu\n", tile.cycles());

    return;
}

void radix2_dit_3(input_window<cint16> *x_in, output_window<cint16> *y_out)
{
    cint16 *x = (cint16 *)x_in->ptr;
    cint16 *y = (cint16 *)y_out->ptr;

    // --------------------------------shuffle--------------------------------

    aie::tile tile = aie::tile::current();
    printf("before shuffle: %llu\n", tile.cycles());

    for (unsigned i = 0; i < N_SWAP_2; i += 2)
    {
        int16 a=swap2[i],b=swap2[i+1];
        cint16 tmp = x[a];
        x[a] = x[b];
        x[b] = tmp;
    }
    for (unsigned i = 0; i < N_SWAP_4; i += 4)
    {
        int16 i0=swap4[i],i1=swap4[i+1],i2=swap4[i+2],i3=swap4[i+3];
        cint16 tmp = x[i3];
        x[i3] = x[i0];
        x[i0] = x[i1];
        x[i1] = x[i2];
        x[i2] = tmp;
    }

    printf("shuffled: %llu\n", tile.cycles());

    // ----------------------------------dit----------------------------------

    for (cint16* p=x;p!=x+N_POINT;p+=MAX_VEC_LEN)
    {
        auto iter=begin_vector<MAX_VEC_LEN>(mat_omg_8);
        auto m=mul(*iter++,*p);

        for (unsigned i=1;i<MAX_VEC_LEN;i++){
            m=mac(m,*iter++,*(p+i));
        }
        store_v(p,m.to_vector<cint16>(MAT_OMG_SHIFT));
    }
    printf("l<=MAX_VEC_LEN: %llu\n", tile.cycles());
    for (unsigned l = MAX_VEC_LEN << 1; l <= N_POINT; l <<= 1)
    {
        butterfly(l, x, x);
        printf("btf l=%u: %llu\n", l, tile.cycles());
    }

    mul_twiddle_factor(x,y_out,(cint16*)tf3);

    printf("dit: %llu\n", tile.cycles());

    return;
}
