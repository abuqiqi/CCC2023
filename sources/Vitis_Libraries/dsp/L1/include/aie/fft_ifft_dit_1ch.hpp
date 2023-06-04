/*
 * Copyright 2022 Xilinx, Inc.
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef _DSPLIB_FFT_IFFT_DIT_1CH_HPP_
#define _DSPLIB_FFT_IFFT_DIT_1CH_HPP_

/*
FFT/IFFT, DIT, single channel.
This file exists to capture the definition of the single channel FFT/iFFT filter kernel class.
The class definition holds defensive checks on parameter range and other legality.
The constructor definition is held in this class because this class must be accessible to graph
level aie compilation.
The main runtime function is captured elsewhere as it contains aie intrinsics which are not
included in aie graph level compilation.
*/

/* Coding conventions
   TT_      template type suffix
   TP_      template parameter suffix
*/

/* Design Notes
*/

#include <adf.h>
#include <vector>

#include "fft_ifft_dit_1ch_traits.hpp"

#ifndef _DSPLIB_FFT_IFFT_DIT_1CH_HPP_DEBUG_
//#define _DSPLIB_FFT_IFFT_DIT_1CH_HPP_DEBUG_
#endif //_DSPLIB_FFT_IFFT_DIT_1CH_HPP_DEBUG_

namespace xf {
namespace dsp {
namespace aie {
namespace fft {
namespace dit_1ch {
//-----------------------------------------------------------------------------------------------------

template <typename TT_DATA = cint16,
          typename TT_OUT_DATA = cint16,
          typename TT_TWIDDLE = cint16,
          typename TT_INTERNAL_DATA = cint16,
          unsigned int TP_POINT_SIZE = 4096,
          unsigned int TP_FFT_NIFFT = 1,
          unsigned int TP_SHIFT = 0,
          unsigned int TP_START_RANK = 0,
          unsigned int TP_END_RANK = 8,
          unsigned int TP_DYN_PT_SIZE = 0,
          unsigned int TP_WINDOW_VSIZE = TP_POINT_SIZE,
          unsigned int TP_ORIG_PAR_POWER = 0>
class stockhamStages {
   private:
   public:
    stockhamStages() {} // Constructor

    static constexpr unsigned int kPointSizePower = fnPointSizePower<TP_POINT_SIZE>();
    static constexpr unsigned int kPointSizePowerCeiled =
        fnCeil<kPointSizePower, 2>(); // round up to 2 to allow loop of r4 stages
    static constexpr unsigned int kPointSizeCeiled = 1 << kPointSizePowerCeiled; // loop is for R4 or R2

    static constexpr unsigned int kOddPower = fnOddPower<TP_POINT_SIZE>();
    static constexpr unsigned int kIntR4Stages =
        (kPointSizePower - 3) / 2; // number of internal radix 4 stags, i.e. not including input or output stages.
    static constexpr unsigned int kIntR2Stages =
        (kPointSizePower -
         2); // The 2 comes from the last 2 ranks being handled differently, being Stockham stage 1 and 2.

    void calc(TT_DATA* __restrict xbuff,
              TT_TWIDDLE** tw_table,
              TT_INTERNAL_DATA* tmp1_buf,
              TT_INTERNAL_DATA* tmp2_buf,
              TT_OUT_DATA* __restrict obuff);
    void calc(TT_DATA* __restrict xbuff,
              TT_TWIDDLE** tw_table,
              TT_INTERNAL_DATA* tmp1_buf,
              TT_INTERNAL_DATA* tmp2_buf,
              TT_OUT_DATA* __restrict obuff,
              int ptSizePwr,
              bool inv);
    void stagePreamble(TT_TWIDDLE** tw_table,
                       TT_INTERNAL_DATA* tmp1_buf,
                       TT_INTERNAL_DATA* tmp2_buf,
                       input_window<TT_DATA>* __restrict inputx,
                       output_window<TT_OUT_DATA>* __restrict outputy);
};

template <typename TT_DATA = cint16,
          typename TT_OUT_DATA = cint16,
          typename TT_TWIDDLE = cint16,
          unsigned int TP_POINT_SIZE = 4096,
          unsigned int TP_FFT_NIFFT = 1,
          unsigned int TP_SHIFT = 0,
          unsigned int TP_START_RANK = 0,
          unsigned int TP_END_RANK = 8,
          unsigned int TP_DYN_PT_SIZE = 0,
          unsigned int TP_WINDOW_VSIZE = TP_POINT_SIZE,
          unsigned int TP_ORIG_PAR_POWER = 0>
class kernelFFTClass {
   private:
   public:
    // Constructor
    kernelFFTClass() {}

    // FFT
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};

//-------------------------------
// Specializations of kernelFFTClass for each point size is required because of the list of twiddles to include
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     4096,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   4096,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     2048,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   2048,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     1024,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   1024,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     512,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   512,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     256,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   256,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     128,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   128,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     64,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   64,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     32,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   32,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};

template <typename TT_DATA,
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<TT_DATA,
                     TT_OUT_DATA,
                     TT_TWIDDLE,
                     16,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef typename std::conditional<std::is_same<TT_DATA, cint16>::value, cint32_t, TT_DATA>::type T_internalDataType;
    stockhamStages<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   T_internalDataType,
                   16,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;
    void kernelFFT(input_window<TT_DATA>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};

// Specialisation for the first kernel with cint16 input. This is the case which requires an explicit extra buffer tmp2
// rather than re-use the input
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     4096,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   4096,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     2048,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   2048,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     1024,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   1024,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     512,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   512,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     256,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   256,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     128,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   128,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     64,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   64,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};
template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     32,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   32,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};

template <typename TT_OUT_DATA,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_ORIG_PAR_POWER>
class kernelFFTClass<cint16,
                     TT_OUT_DATA,
                     cint16,
                     16,
                     TP_FFT_NIFFT,
                     TP_SHIFT,
                     TP_START_RANK,
                     TP_END_RANK,
                     TP_DYN_PT_SIZE,
                     TP_WINDOW_VSIZE,
                     TP_ORIG_PAR_POWER> {
   public:
    typedef cint32_t T_internalDataType;
    stockhamStages<cint16,
                   TT_OUT_DATA,
                   cint16,
                   T_internalDataType,
                   16,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        stages;

    void kernelFFT(input_window<cint16>* __restrict inputx, output_window<TT_OUT_DATA>* __restrict outputy);
};

//-----------------------------------------------------------------------------------------------------
// Top level single kernel specialization.
template <typename TT_DATA, // input
          typename TT_OUT_DATA,
          typename TT_TWIDDLE,
          unsigned int TP_POINT_SIZE,
          unsigned int TP_FFT_NIFFT,
          unsigned int TP_SHIFT,
          unsigned int TP_START_RANK,
          unsigned int TP_END_RANK,
          unsigned int TP_DYN_PT_SIZE,
          unsigned int TP_WINDOW_VSIZE = TP_POINT_SIZE,
          unsigned int TP_ORIG_PAR_POWER = 0>
class fft_ifft_dit_1ch {
   private:
    // Parameter value defensive and legality checks
    static_assert(fnCheckDataType<TT_DATA>(), "ERROR: TT_IN_DATA is not a supported type");
    static_assert(fnCheckDataType<TT_OUT_DATA>(), "ERROR: TT_OUT_DATA is not a supported type");
    static_assert(fnCheckDataIOType<TT_DATA, TT_OUT_DATA>(), "ERROR: TT_OUT_DATA is not a supported type");
    static_assert(fnCheckTwiddleType<TT_TWIDDLE>(), "ERROR: TT_TWIDDLE is not a supported type");
    static_assert(fnCheckDataTwiddleType<TT_DATA, TT_TWIDDLE>(), "ERROR: TT_TWIDDLE is incompatible with data type");
    static_assert(fnCheckPointSize<TP_POINT_SIZE>(),
                  "ERROR: TP_POINT_SIZE is not a supported value {16, 32, 64, ..., 4096}");
    static_assert(TP_FFT_NIFFT == 0 || TP_FFT_NIFFT == 1, "ERROR: TP_FFT_NIFFT must be 0 (reverse) or 1 (forward)");
    static_assert(fnCheckShift<TP_SHIFT>(), "ERROR: TP_SHIFT is out of range (0 to 60)");
    static_assert(fnCheckShiftFloat<TT_DATA, TP_SHIFT>(),
                  "ERROR: TP_SHIFT is ignored for data type cfloat so must be set to 0");
    static_assert(TP_WINDOW_VSIZE % TP_POINT_SIZE == 0, "ERROR: TP_WINDOW_SIZE must be a multiple of TP_POINT_SIZE");
    static_assert(TP_WINDOW_VSIZE / TP_POINT_SIZE >= 1, "ERROR: TP_WINDOW_SIZE must be a multiple of TP_POINT_SIZE");
    static_assert((TP_DYN_PT_SIZE == 0) || (TP_POINT_SIZE != 16),
                  "ERROR: Dynamic point size is not supported for only one legal size (16)");

    kernelFFTClass<TT_DATA,
                   TT_OUT_DATA,
                   TT_TWIDDLE,
                   TP_POINT_SIZE,
                   TP_FFT_NIFFT,
                   TP_SHIFT,
                   TP_START_RANK,
                   TP_END_RANK,
                   TP_DYN_PT_SIZE,
                   TP_WINDOW_VSIZE,
                   TP_ORIG_PAR_POWER>
        m_fftKernel; // Kernel for FFT

   public:
    // Constructor
    fft_ifft_dit_1ch() {}

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(fft_ifft_dit_1ch::fftMain); }
    // FFT
    void fftMain(input_window<TT_DATA>* __restrict inWindow, output_window<TT_OUT_DATA>* __restrict outWindow);
};
}
}
}
}
}
#endif // _DSPLIB_FFT_IFFT_DIT_1CH_HPP_
