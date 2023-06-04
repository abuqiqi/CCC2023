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
/*
This file holds the body of the kernel class for the Asymmetric Interpolation FIR.
Unlike single rate implementations, this interpolation FIR calculates sets of output
vectors in parallel such that the number of lanes in total is a multiple of the
interpolation factor.

Coding conventions
  TT_      template type suffix
  TP_      template parameter suffix
*/

#pragma once
#include <adf.h>

#define __AIE_API_USE_NATIVE_1024B_VECTOR__
#include "aie_api/aie_adf.hpp"

#include "kernel_api_utils.hpp"

#include "fir_interpolate_asym.hpp"
#include "fir_interpolate_asym_utils.hpp"
#ifndef Y_BUFFER
#define Y_BUFFER yd
#endif

#ifndef ACC_BUFFER
#define ACC_BUFFER bm0
#endif

#ifndef OUT_BUFFER
#define OUT_BUFFER wr0
#endif

#ifndef OUT_BUFFER_SHORT
#define OUT_BUFFER_SHORT vrh0
#endif

namespace xf {
namespace dsp {
namespace aie {
namespace fir {
namespace interpolate_asym {
// What follows are templatized functions to precompute values used by the run-time filter
// function so as to minimise run-time computation and instead look-up precomputed values.
// Specifically the xoffsets arguments to the MAC intrinsic central to this algorithm are calculated here.
template <unsigned TP_INTERPOLATE_FACTOR, typename TT_DATA, typename TT_COEFF, int TP_MODIFY_COEFF_OFFSET = 0>
INLINE_DECL constexpr std::array<int32, TP_INTERPOLATE_FACTOR> fnInitStarts() {
    constexpr unsigned int m_kXbuffSize = 128; // kXbuffSize in Bytes (1024bit) - const for all data/coeff types
    constexpr unsigned int m_kDataRegVsize = m_kXbuffSize / sizeof(TT_DATA); // sbuff size in Bytes
    constexpr unsigned int m_kLanes = fnNumLanesIntAsym<TT_DATA, TT_COEFF>();
    std::array<int32, TP_INTERPOLATE_FACTOR> ret = {};
    for (int phase = 0; phase < TP_INTERPOLATE_FACTOR; ++phase) {
        int effPhase = (phase * m_kLanes + TP_MODIFY_COEFF_OFFSET); // the coefficient phase for the given lane changes
                                                                    // by TP_MODIFY_COEFF_OFFSET due to SSR.
        int coeffPhase = FLOOR(effPhase, TP_INTERPOLATE_FACTOR);
        ret[phase] =
            coeffPhase > 0
                ? (coeffPhase % m_kDataRegVsize)
                : -1 * (-1 * coeffPhase % m_kDataRegVsize); // expecting ret[phase] to have same sign as coeffPhase
    }
    return ret;
};

template <unsigned TP_INTERPOLATE_FACTOR, typename TT_DATA, typename TT_COEFF, int TP_MODIFY_COEFF_OFFSET = 0>
INLINE_DECL constexpr std::array<uint64, TP_INTERPOLATE_FACTOR> fnInitOffsets() {
    constexpr unsigned int m_kLanes =
        fnNumLanesIntAsym<TT_DATA, TT_COEFF>(); // kMaxMacs/(sizeof(TT_DATA)*sizeof(TT_COEFF)*m_kColumns); //number of
                                                // operations in parallel of this type combinations that the vector
                                                // processor can do.
    constexpr unsigned int bitsInNibble = 4;
    constexpr unsigned int m_kXbuffSize = 128; // kXbuffSize in Bytes (1024bit) - const for all data/coeff types
    constexpr unsigned int m_kDataRegVsize = m_kXbuffSize / sizeof(TT_DATA); // sbuff size in Bytes
    std::array<uint64, TP_INTERPOLATE_FACTOR> ret = {};
    uint64 dataEntry = 0;
    for (int phase = 0; phase < TP_INTERPOLATE_FACTOR; ++phase) {
        dataEntry = 0;
        int effPhase = (phase * m_kLanes + TP_MODIFY_COEFF_OFFSET); // the coefficient phase for the given lane changes
                                                                    // by TP_MODIFY_COEFF_OFFSET due to SSR.
        int coeffPhase = FLOOR(effPhase, TP_INTERPOLATE_FACTOR);
        int initStartForPhase =
            coeffPhase > 0
                ? (coeffPhase % m_kDataRegVsize)
                : -1 * (-1 * coeffPhase % m_kDataRegVsize); // expecting ret[phase] to have same sign as coeffPhase
        for (int lane = 0; lane < m_kLanes; ++lane) {
            effPhase = (phase * m_kLanes + lane + TP_MODIFY_COEFF_OFFSET);
            coeffPhase = FLOOR(effPhase, TP_INTERPOLATE_FACTOR);
            int initOffset = coeffPhase > 0 ? (coeffPhase % m_kDataRegVsize) : -1 * (-1 * coeffPhase % m_kDataRegVsize);
            dataEntry += (initOffset - initStartForPhase) << (lane * bitsInNibble);
        }
        // Note that m_dataOffsets does not vary by op (m_dataStarts does), but code is left this way to show derivation
        // of values
        ret[phase] = dataEntry;
    }
    return ret;
};

// FIR function
//-----------------------------------------------------------------------------------------------------
//#TEMPLATE_FUNCTION_DEFINITION
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filterKernel(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                     T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    windowBroadcast<TT_DATA, TP_INPUT_WINDOW_VSIZE + fnFirMargin<TP_FIR_LEN / TP_INTERPOLATE_FACTOR, TT_DATA>(),
                    TP_API>(inInterface, outInterface);
    filterSelectArch(inInterface, outInterface);
}
// FIR function
//-----------------------------------------------------------------------------------------------------
//#TEMPLATE_FUNCTION_DEFINITION
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filterKernel(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                     T_outputIF<TP_CASC_OUT, TT_DATA> outInterface,
                                                     const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    windowBroadcast<TT_DATA, TP_INPUT_WINDOW_VSIZE + fnFirMargin<TP_FIR_LEN / TP_INTERPOLATE_FACTOR, TT_DATA>(),
                    TP_API>(inInterface, outInterface);
    isUpdateRequired = rtpCompare<TT_COEFF, TP_COEFF_PHASES_LEN>(inTaps, m_rawInTaps);

    sendRtpTrigger(isUpdateRequired, outInterface);
    if (isUpdateRequired) { // Coefficients have changed
        bufferReload<TT_DATA, TT_COEFF, TP_COEFF_PHASES_LEN>(inTaps, m_rawInTaps, outInterface);
        firReload<TP_COEFF_PHASE, TP_COEFF_PHASE_OFFSET, TP_COEFF_PHASES, TP_COEFF_PHASES_LEN>(inTaps);
    }
    filterSelectArch(inInterface, outInterface);
}

// Asymmetric Fractional Interpolation FIR Kernel Function - overloaded (not specialised)
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filterKernelRtp(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                        T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    windowBroadcast<TT_DATA, TP_INPUT_WINDOW_VSIZE + fnFirMargin<TP_FIR_LEN / TP_INTERPOLATE_FACTOR, TT_DATA>(),
                    TP_API>(inInterface, outInterface);
    isUpdateRequired = getRtpTrigger(); // 0 - equal, 1 - not equal

    sendRtpTrigger(isUpdateRequired, outInterface);
    if (isUpdateRequired) { // Coefficients have changed
        bufferReload<TT_DATA, TT_COEFF, TP_COEFF_PHASES_LEN>(inInterface, m_rawInTaps, outInterface);
        firReload<TP_COEFF_PHASE, TP_COEFF_PHASE_OFFSET, TP_COEFF_PHASES, TP_COEFF_PHASES_LEN>(m_rawInTaps);
    }
    filterSelectArch(inInterface, outInterface);
}

//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filterSelectArch(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                         T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    set_rnd(TP_RND);
    set_sat();

    windowReset<TT_DATA, TP_CASC_IN, TP_DUAL_IP, TP_API>(inInterface);
    // This function hides exposure of the implementation choice from the user.
    if
        constexpr(TP_API == USE_STREAM_API) { filterStream(inInterface, outInterface); }
    else if
        constexpr(m_kArch == kArchPhaseSeries) { filter_impl1(inInterface, outInterface); }
    else if
        constexpr(m_kArch == kArchPhaseParallel) { filter_impl2(inInterface, outInterface); }
    else {
        filterIncr(inInterface, outInterface);
    }
    windowRelease(inInterface);
}

// Implementation 1, Here, each of the phases is calculated in series to avoid pulling and pushing
// the accumulator to the stack.
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filterStream(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                     T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    set_rnd(TP_RND);
    set_sat();
    T_buff_256b<TT_COEFF>* __restrict coeff2 = (T_buff_256b<TT_COEFF>*)m_internalTaps;
    TT_COEFF* coeff = (TT_COEFF*)m_internalTaps;

    T_buff_256b<TT_COEFF> coe0;                                       // register for coeff values.
    T_buff_1024b<TT_DATA>* ptr_delay = (T_buff_1024b<TT_DATA>*)delay; // heap storage pointer

    T_buff_1024b<TT_DATA> sbuff = *ptr_delay; // cache register for data values entering MAC
    using buf_type = typename T_buff_1024b<TT_DATA>::v_type;
    buf_type chess_storage(Y_BUFFER) sbuffTmp;

    using outTypeInt = typename T_outValIntAsym<TT_DATA, TT_COEFF>::v_type chess_storage(OUT_BUFFER);
    using outTypeIntShort = typename T_outValIntAsym<TT_DATA, TT_COEFF>::v_type chess_storage(OUT_BUFFER_SHORT);
    using outType =
        typename std::conditional<(std::is_same<TT_DATA, cint32>::value) && (std::is_same<TT_COEFF, cint32>::value),
                                  outTypeInt, outTypeIntShort>::type;

    outType outValTmp;

    T_accIntAsym<TT_DATA, TT_COEFF> acc, debugAcc;
    T_outValIntAsym<TT_DATA, TT_COEFF> outVal;
    T_buff_256b<TT_DATA> readData;

    unsigned int index = startIndex;
    int xstart;
    unsigned int LCMPhase = 0;
    unsigned int dataLoaded, dataNeeded; // In terms of register locations, not data samples
    constexpr static std::array<int32, TP_INTERPOLATE_FACTOR> m_dataStarts =
        fnInitStarts<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF, TP_MODIFY_MARGIN_OFFSET>();
    constexpr static std::array<uint64, TP_INTERPOLATE_FACTOR> m_dataOffsets =
        fnInitOffsets<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF, TP_MODIFY_MARGIN_OFFSET>();
    int loopSize = (m_kLsize / m_kRepeatFactor);
    TT_COEFF* m_internalTapsCopy = (TT_COEFF*)m_internalTaps; // points to m_internalTaps[0][0][0][0]

    if
        constexpr(TP_CASC_LEN > 1) {
            if (doInit == 1) {
                dataLoaded = 0;
                LCMPhase = phaseOffset;
                index = 2;
                for (unsigned i = 0; i < streamInitNullAccs; i++)
                    chess_prepare_for_pipelining chess_loop_range(streamInitNullAccs, ) {
                        acc = readCascade(inInterface, acc);

                        writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);
                    }
#pragma unroll(GUARD_ZERO(CEIL(streamInitNullAccs, m_kPhases) - streamInitNullAccs))
                for (unsigned i = streamInitNullAccs; i < CEIL(streamInitNullAccs, m_kPhases); i++) {
                    m_internalTapsCopy = chess_copy(m_internalTapsCopy);
                    coeff = m_internalTapsCopy + (LCMPhase * m_kLanes * m_kColumns * m_kNumOps);
                    coe0 = *(T_buff_256b<TT_COEFF>*)coeff;
                    dataNeeded = ((int)(((m_kLanes * ((i % TP_INTERPOLATE_FACTOR) + 1)) - 1) / TP_INTERPOLATE_FACTOR) +
                                  1 + streamInitNullAccs / m_kPhases * totalNeededDataPerStrobe) -
                                 initPhaseDataAccesses;

                    if (dataNeeded > dataLoaded) {
                        if
                            constexpr(TP_DUAL_IP == DUAL_IP_SINGLE) {
                                readStream128(sbuff, index, inInterface);
                                sbuffTmp = sbuff.val;
                                sbuff.val = sbuffTmp;
                                index++;
                            }
                        else {
                            readStream256(sbuff, index, inInterface);
                            index++;
                        }
                        dataLoaded += sizeOf1Read;
                    }
                    xstart = m_kXStart + streamInitNullAccs / m_kPhases * totalNeededDataPerStrobe +
                             m_dataStarts[i % (TP_INTERPOLATE_FACTOR)];
                    acc = readCascade(inInterface, acc);

                    // accTmp = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                    // m_dataOffsets[LCMPhase], xstart).val;
                    // acc.val = accTmp;
                    acc = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                         m_dataOffsets[LCMPhase], xstart);

#pragma unroll(GUARD_ZERO((m_kNumOps - 1)))
                    for (int op = 1; op < m_kNumOps; ++op) {
                        xstart += m_kColumns;
                        coeff = m_internalTapsCopy + (LCMPhase * m_kLanes * m_kColumns * m_kNumOps) +
                                (op * m_kLanes * m_kColumns);
                        coe0 = *(T_buff_256b<TT_COEFF>*)coeff;

                        acc = macIntAsym(acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes, m_dataOffsets[LCMPhase],
                                         xstart);
                    }

                    outValTmp = shiftAndSaturateIntAsym(acc, TP_SHIFT).val;
                    outVal.val = outValTmp;

                    // Write cascade. Do nothing if cascade not present.
                    writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);

                    writeOutput<TT_DATA, TT_COEFF, TP_NUM_OUTPUTS, TP_API>(outInterface, outVal, i % 2);
                    LCMPhase = (LCMPhase + 1);
                    if (LCMPhase == m_kLCMPhases) { // a mod function without the division.TODO optimize TP_INTERPOLATE
                                                    // to prime factor
                        LCMPhase = 0;
                    }
                }

#pragma unroll(GUARD_ZERO(streamInitAccs))
                for (unsigned strobe = 0; strobe < (streamInitAccs); strobe++) {
#pragma unroll(m_kPhases)
                    // The phase loop effectively multiplies the number of lanes in use to ensures that
                    // an integer number of interpolation polyphases are calculated
                    for (int phase = 0; phase < m_kPhases; ++phase) {
                        m_internalTapsCopy = chess_copy(m_internalTapsCopy);
                        coeff = (m_internalTapsCopy + (LCMPhase * m_kLanes * m_kColumns * m_kNumOps));
                        coe0 = *(T_buff_256b<TT_COEFF>*)coeff;
                        // coe0 = null_buff_1024b<TT_DATA>();
                        dataNeeded = ((strobe + streamInitNullStrobes) * totalNeededDataPerStrobe +
                                      (int)((m_kLanes * (phase + 1) - 1) / TP_INTERPOLATE_FACTOR) + 1) -
                                     initPhaseDataAccesses;

                        if (dataNeeded > dataLoaded) {
                            if
                                constexpr(TP_DUAL_IP == DUAL_IP_SINGLE) {
                                    readStream128(sbuff, index, inInterface);
                                    sbuffTmp = sbuff.val;
                                    sbuff.val = sbuffTmp;
                                    index++;
                                }
                            else {
                                readStream256(sbuff, index, inInterface);
                                index++;
                            }
                            dataLoaded += sizeOf1Read;
                        }

                        xstart = m_kXStart + (strobe + streamInitNullStrobes) * m_kLanes + m_dataStarts[phase];
                        acc = readCascade(inInterface, acc);

                        acc = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                             m_dataOffsets[LCMPhase], xstart);

#pragma unroll(GUARD_ZERO((m_kNumOps - 1)))
                        for (int op = 1; op < m_kNumOps; ++op) {
                            xstart += m_kColumns;
                            // coeff =  ((T_buff_256b<TT_COEFF> *)&m_internalTaps[LCMPhase][op][0][0]);
                            coeff = m_internalTapsCopy + (LCMPhase * m_kLanes * m_kColumns * m_kNumOps) +
                                    (op * m_kLanes * m_kColumns);
                            coe0 = *(T_buff_256b<TT_COEFF>*)coeff;

                            acc = macIntAsym(acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes, m_dataOffsets[LCMPhase],
                                             xstart);
                        }

                        outValTmp = shiftAndSaturateIntAsym(acc, TP_SHIFT).val;
                        outVal.val = outValTmp;

                        // Write cascade. Do nothing if cascade not present.
                        writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);

                        writeOutput<TT_DATA, TT_COEFF, TP_NUM_OUTPUTS, TP_API>(outInterface, outVal, phase % 2);
                        LCMPhase = (LCMPhase + 1);
                        if (LCMPhase == m_kLCMPhases) { // a mod function without the division.TODO optimize
                                                        // TP_INTERPOLATE to prime factor
                            LCMPhase = 0;
                        }
                    }
                }
                loopSize -= CEIL(streamInitNullAccs, m_kRepFactPhases) / (m_kRepeatFactor * m_kPhases);
            }
        }
    doInit = 0;

    coeff = m_internalTapsCopy;
    // loop through window, computing a vector of output for each iteration.
    for (unsigned i = 0; i < loopSize; i++)
        chess_prepare_for_pipelining chess_pipeline_non_leaf_loop_solution(4) chess_loop_range(
            (m_kLsize / m_kRepeatFactor) - (CEIL(streamInitNullAccs, m_kRepFactPhases) / (m_kRepeatFactor * m_kPhases)),
            (m_kLsize / m_kRepeatFactor)) {
            dataLoaded = 0;
            LCMPhase = 0;
            if
                constexpr(streamInitNullAccs == 0) { index = startIndex; }
            else {
                index = startIndexCasc;
            }
#pragma unroll(m_kRepeatFactor)
            for (int strobe = 0; strobe < m_kRepeatFactor; strobe++) {
#pragma unroll(m_kPhases)
                // The phase loop effectively multiplies the number of lanes in use to ensures that
                // an integer number of interpolation polyphases are calculated
                for (int phase = 0; phase < m_kPhases; ++phase) {
                    coeff = chess_copy(coeff);
                    coe0 = *(T_buff_256b<TT_COEFF>*)coeff;
                    coeff = coeff + (m_kLanes * m_kColumns);

                    if
                        constexpr(streamInitNullAccs == 0) {
                            dataNeeded = strobe * ((int)((m_kLanes * (m_kPhases)-1) / TP_INTERPOLATE_FACTOR) + 1) +
                                         (int)((m_kLanes * (phase + 1) - 1) / TP_INTERPOLATE_FACTOR) + 1;
                        }
                    else {
                        dataNeeded = strobe * ((int)((m_kLanes * (m_kPhases)-1) / TP_INTERPOLATE_FACTOR) + 1) +
                                     (int)((m_kLanes * (phase + 1) - 1) / TP_INTERPOLATE_FACTOR) + 1;
                    }

                    if (dataNeeded > dataLoaded) {
                        if
                            constexpr(TP_DUAL_IP == DUAL_IP_SINGLE) {
                                readStream128(sbuff, index, inInterface);
                                sbuffTmp = sbuff.val;
                                sbuff.val = sbuffTmp;
                                index++;
                            }
                        else {
                            readStream256(sbuff, index, inInterface);
                            index++;
                        }
                        dataLoaded += sizeOf1Read;
                    }

                    xstart = m_kXStart + strobe * m_kLanes + m_dataStarts[phase];
                    acc = readCascade(inInterface, acc);

                    acc = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                         m_dataOffsets[LCMPhase], xstart);

#pragma unroll(GUARD_ZERO((m_kNumOps - 1)))
                    for (int op = 1; op < m_kNumOps; ++op) {
                        xstart += m_kColumns;
                        coe0 = *(T_buff_256b<TT_COEFF>*)coeff;
                        if
                            constexpr(m_kPhases == m_kLCMPhases) {
                                if (op == m_kNumOps - 1) {
                                    if (phase == m_kPhases - 1) {
                                        coeff = coeff -
                                                (m_kPhases * m_kLanes * m_kColumns * m_kNumOps - m_kLanes * m_kColumns);
                                    } else {
                                        coeff = coeff + (m_kLanes * m_kColumns);
                                    }
                                } else {
                                    coeff = coeff + (m_kLanes * m_kColumns);
                                }
                            }
                        else {
                            if (op == m_kNumOps - 1) {
                                if (LCMPhase == m_kLCMPhases - 1) {
                                    coeff = coeff -
                                            (m_kLCMPhases * m_kLanes * m_kColumns * m_kNumOps - m_kLanes * m_kColumns);
                                } else {
                                    coeff = coeff + (m_kLanes * m_kColumns);
                                }
                            } else {
                                coeff = coeff + (m_kLanes * m_kColumns);
                            }
                        }
                        acc = macIntAsym(acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes, m_dataOffsets[LCMPhase],
                                         xstart);
                    }
                    if
                        constexpr(m_kNumOps == 1) {
                            if (LCMPhase == m_kLCMPhases - 1) {
                                coeff = coeff - (m_kLCMPhases * m_kLanes * m_kColumns * m_kNumOps);
                            }
                        }

                    // Write cascade. Do nothing if cascade not present.
                    writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);

                    outValTmp = shiftAndSaturateIntAsym(acc, TP_SHIFT).val;
                    outVal.val = outValTmp;

                    writeOutput<TT_DATA, TT_COEFF, TP_NUM_OUTPUTS, TP_API>(outInterface, outVal, phase % 2);
                    LCMPhase = (LCMPhase + 1);
                    if (LCMPhase == m_kLCMPhases) { // a mod function without the division.TODO optimize TP_INTERPOLATE
                                                    // to prime factor
                        LCMPhase = 0;
                    }
                }
            }
        }
    doInit = 0;
    *ptr_delay = sbuff;
}

// Implementation 1, Here, each of the phases is calculated in series to avoid pulling and pushing
// the accumulator to the stack.
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filter_impl1(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                     T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    set_rnd(TP_RND);
    set_sat();

    T_buff_256b<TT_COEFF>* __restrict coeff = (T_buff_256b<TT_COEFF>*)m_internalTaps;
    TT_COEFF* m_internalTapsCopy = (TT_COEFF*)m_internalTaps; // points to m_internalTaps[0][0][0][0]
    T_buff_256b<TT_COEFF> coe0;                               // register for coeff values.
    T_buff_1024b<TT_DATA> sbuff;                              // cache register for data values entering MAC
    T_accIntAsym<TT_DATA, TT_COEFF> acc, debugAcc;
    T_outValIntAsym<TT_DATA, TT_COEFF> outVal;
    T_buff_256b<TT_DATA> readData;
    unsigned int dataLoaded, dataNeeded; // In terms of register locations, not data samples
    unsigned int numDataLoads;
    unsigned int LCMPhase = 0;
    unsigned int xstart = 0;
    uint64 xoffsets = 0;
    // m_dataOffsets allows offsets to be precalculated in the constructor then simply looked up according to loop
    // indices during runtime.
    constexpr static std::array<int32, TP_INTERPOLATE_FACTOR> m_dataStarts =
        fnInitStarts<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF>();
    constexpr static std::array<uint64, TP_INTERPOLATE_FACTOR> m_dataOffsets =
        fnInitOffsets<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF>();

    window_incr(inInterface.inWindow, m_kDataWindowOffset); // move input data pointer past the margin padding
    // loop through window, computing a vector of output for each iteration.
    for (unsigned i = 0; i < m_kLsize; i++) chess_prepare_for_pipelining chess_loop_range(m_kLsize, ) {
#pragma unroll(m_kPhases)
            // The phase loop effectively multiplies the number of lanes in use to ensures that
            // an integer number of interpolation polyphases are calculated
            for (int phase = 0; phase < m_kPhases; ++phase) {
                m_internalTapsCopy = chess_copy(m_internalTapsCopy);
                coeff = (T_buff_256b<TT_COEFF>*)(m_internalTapsCopy + (LCMPhase * m_kLanes * m_kColumns * m_kNumOps));
                coe0 = *coeff;

                numDataLoads = 0;
                dataLoaded = 0;

                // preamble, load data from window into register
                dataNeeded = m_kDataBuffXOffset * TP_INTERPOLATE_FACTOR + m_kLanes * (phase + 1) +
                             (m_kColumns - 1) * TP_INTERPOLATE_FACTOR;

// load the data registers with enough data for the initial MAC(MUL)
#pragma unroll(m_kInitialLoads)
                for (int initLoads = 0; initLoads < m_kInitialLoads; ++initLoads) {
                    readData = window_readincr_256b<TT_DATA>(inInterface.inWindow);
                    sbuff.val = upd_w(sbuff.val, numDataLoads, readData.val);
                    dataLoaded += m_kDataLoadVsize * TP_INTERPOLATE_FACTOR; // in effect, since data is duplicated
                    numDataLoads++;
                }
                xstart = m_kDataBuffXOffset + m_dataStarts[phase];
                acc = readCascade<TT_DATA, TT_COEFF>(inInterface, acc);
                // Perform the first term (or 2 or 4) of the FIR polynomial.
                acc = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                     m_dataOffsets[LCMPhase], xstart);

// loop to work through the operations in the FIR polynomial
#pragma unroll(GUARD_ZERO((m_kNumOps - 1)))
                for (int op = 1; op < m_kNumOps; ++op) {
                    dataNeeded += m_kColumns * TP_INTERPOLATE_FACTOR;
                    if (dataNeeded > dataLoaded) {
                        readData = window_readincr_256b<TT_DATA>(inInterface.inWindow);
                        sbuff.val =
                            upd_w(sbuff.val, numDataLoads % m_kDataLoadsInReg, readData.val); // ____|____|XX..|XX...
                        dataLoaded += m_kDataLoadVsize * TP_INTERPOLATE_FACTOR;
                        numDataLoads++;
                    }
                    // coeff =  ((T_buff_256b<TT_COEFF> *)&m_internalTaps[LCMPhase][op][0][0]);
                    coeff =
                        (T_buff_256b<TT_COEFF>*)(m_internalTapsCopy + (LCMPhase * m_kLanes * m_kColumns * m_kNumOps) +
                                                 (op * m_kLanes * m_kColumns));
                    coe0 = *coeff;
                    xstart += m_kColumns;
                    // perform the MAC, i.e. cacculate some terms of the FIR polynomial
                    acc =
                        macIntAsym(acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes, m_dataOffsets[LCMPhase], xstart);
                }

                // Write cascade. Do nothing if cascade not present.
                writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);

                // convert the accumulator type into an integer type for output, downshift and apply any rounding and
                // saturation.
                outVal = shiftAndSaturateIntAsym(acc, TP_SHIFT);

                // Write to output window
                writeWindow<TT_DATA, TT_COEFF, TP_NUM_OUTPUTS>(outInterface, outVal);

                LCMPhase = (LCMPhase + 1);
                if (LCMPhase ==
                    m_kLCMPhases) { // a mod function without the division.TODO optimize TP_INTERPOLATE to prime factor
                    LCMPhase = 0;
                }

                // take data pointer back to next start point.
                window_decr(inInterface.inWindow,
                            (m_kDataLoadVsize * numDataLoads)); // return read pointer to start of next chunk of window.
            }
            window_incr(inInterface.inWindow, m_kLanes); // after all phases, one m_kLanes of input will be consumed.
        }
};

// Implementation 2 is used where there are enough accumulators for the phases, as this implementation requires
// fewer data loads overall.
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filter_impl2(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                     T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    set_rnd(TP_RND);
    set_sat();

    T_buff_256b<TT_COEFF>* __restrict coeff = (T_buff_256b<TT_COEFF>*)m_internalTaps;
    T_buff_256b<TT_COEFF> coe0;  // register for coeff values.
    T_buff_1024b<TT_DATA> sbuff; // undef_buff_1024b<TT_DATA>();
    T_accIntAsym<TT_DATA, TT_COEFF> acc[TP_INTERPOLATE_FACTOR], debugAcc;
    T_outValIntAsym<TT_DATA, TT_COEFF> outVal;
    T_buff_256b<TT_DATA> readData;
    // Note on dataLoaded, dataNeeded. To avoid fractions, both are multiplied by TP_INTERPOLATE_FACTOR.
    // eg. with 8 lanes and a factor of 3, the first operation requires 2.67 samples. Instead, this is accounted
    // as 8 needed, so if 3 samples were loaded, this would be accounted as dataLoaded = 9.
    unsigned int dataLoaded, dataNeeded; // In terms of register locations, not data samples
    unsigned int numDataLoads;
    unsigned int LCMPhase = 0;
    unsigned int xstarts[m_kPhases];

    // m_dataOffsets allows offsets to be precalculated in the constructor then simply looked up according to loop
    // indices during runtime.
    constexpr static std::array<int32, TP_INTERPOLATE_FACTOR> m_dataStarts =
        fnInitStarts<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF>();
    constexpr static std::array<uint64, TP_INTERPOLATE_FACTOR> m_dataOffsets =
        fnInitOffsets<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF>();

    window_incr(inInterface.inWindow, m_kDataWindowOffset); // move input data pointer past the margin padding

    // loop through window, computing a vector of output for each iteration.
    for (unsigned i = 0; i < m_kLsize; ++i) chess_prepare_for_pipelining chess_loop_range(m_kLsize, ) {
            numDataLoads = 0;
            dataLoaded = 0;
#pragma unroll(m_kInitialLoads)
            for (int initLoads = 0; initLoads < m_kInitialLoads; ++initLoads) {
                readData = window_readincr_256b<TT_DATA>(inInterface.inWindow);
                sbuff.val = upd_w(sbuff.val, numDataLoads, readData.val); // 00++|____|____|____
                dataLoaded += m_kDataLoadVsize * TP_INTERPOLATE_FACTOR;   // in effect, since data is duplicated
                numDataLoads++;
            }

            // preamble, load data from window into register
            dataNeeded = m_kDataBuffXOffset * TP_INTERPOLATE_FACTOR + m_kLanes * m_kPhases + m_kColumns - 1;

            LCMPhase = 0;
#pragma unroll(m_kPhases)
            for (int phase = 0; phase < m_kPhases; ++phase) {
                coeff = ((T_buff_256b<TT_COEFF>*)&m_internalTaps[LCMPhase][0][0][0]);
                coe0 = *coeff;
                xstarts[phase] = m_kDataBuffXOffset + m_dataStarts[phase];

                acc[phase] = readCascade<TT_DATA, TT_COEFF>(inInterface, acc);
                acc[phase] = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                            m_dataOffsets[LCMPhase], xstarts[phase]);
                LCMPhase = (LCMPhase + 1);
                if (LCMPhase ==
                    m_kLCMPhases) { // a mod function without the division.TODO optimize TP_INTERPOLATE to prime factor
                    LCMPhase = 0;
                }
            }
#pragma unroll(GUARD_ZERO((m_kNumOps - 1)))
            for (int op = 1; op < m_kNumOps; ++op) {
                dataNeeded += m_kColumns * TP_INTERPOLATE_FACTOR;
                if (dataNeeded > dataLoaded) {
                    readData = window_readincr_256b<TT_DATA>(inInterface.inWindow);
                    sbuff.val =
                        upd_w(sbuff.val, numDataLoads % m_kDataLoadsInReg, readData.val); // ____|____|XX..|XX...
                    dataLoaded += m_kDataLoadVsize * TP_INTERPOLATE_FACTOR;
                    numDataLoads++;
                }
                LCMPhase = 0;
                for (int phase = 0; phase < m_kPhases; ++phase) {
                    xstarts[phase] = (xstarts[phase] + m_kColumns); // modulo count
                    coeff = ((T_buff_256b<TT_COEFF>*)&m_internalTaps[LCMPhase][op][0][0]);
                    coe0 = *coeff;
                    acc[phase] = macIntAsym(acc[phase], sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                            m_dataOffsets[LCMPhase], xstarts[phase]);

                    LCMPhase = (LCMPhase + 1);
                    if (LCMPhase == m_kLCMPhases) { // a mod function without the division.TODO optimize TP_INTERPOLATE
                                                    // to prime factor
                        LCMPhase = 0;
                    }
                }
            }
            for (int phase = 0; phase < m_kPhases; ++phase) {
                // Write cascade. Do nothing if cascade not present.
                writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);
                outVal = shiftAndSaturateIntAsym(acc[phase], TP_SHIFT);
                // Write to output window
                writeWindow<TT_DATA, TT_COEFF, TP_NUM_OUTPUTS>(outInterface, outVal);
            }

            // take data pointer back to next start point.
            window_decr(inInterface.inWindow, (m_kDataLoadVsize * numDataLoads) -
                                                  m_kLanes); // return read pointer to start of next chunk of window.
        }
};

// Implementation 1, Here, each of the phases is calculated in series to avoid pulling and pushing
// the accumulator to the stack.
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
INLINE_DECL void
kernelFilterClass<TT_DATA,
                  TT_COEFF,
                  TP_FIR_LEN,
                  TP_INTERPOLATE_FACTOR,
                  TP_SHIFT,
                  TP_RND,
                  TP_INPUT_WINDOW_VSIZE,
                  TP_CASC_IN,
                  TP_CASC_OUT,
                  TP_FIR_RANGE_LEN,
                  TP_KERNEL_POSITION,
                  TP_CASC_LEN,
                  TP_USE_COEFF_RELOAD,
                  TP_DUAL_IP,
                  TP_NUM_OUTPUTS,
                  TP_API,
                  TP_MODIFY_MARGIN_OFFSET,
                  TP_COEFF_PHASE,
                  TP_COEFF_PHASE_OFFSET,
                  TP_COEFF_PHASES,
                  TP_COEFF_PHASES_LEN>::filterIncr(T_inputIF<TP_CASC_IN, TT_DATA, TP_DUAL_IP> inInterface,
                                                   T_outputIF<TP_CASC_OUT, TT_DATA> outInterface) {
    set_rnd(TP_RND);
    set_sat();

    T_buff_256b<TT_COEFF>* __restrict coeff = (T_buff_256b<TT_COEFF>*)m_internalTaps;
    T_buff_256b<TT_COEFF> coe0;  // register for coeff values.
    T_buff_1024b<TT_DATA> sbuff; // cache register for data values entering MAC
    T_accIntAsym<TT_DATA, TT_COEFF> acc, debugAcc;
    T_outValIntAsym<TT_DATA, TT_COEFF> outVal;
    T_buff_256b<TT_DATA> readData;
    unsigned int dataLoaded, dataNeeded; // In terms of register locations, not data samples
    unsigned int numDataLoads;
    unsigned int LCMPhase = 0;
    unsigned int xstart = 0;
    uint64 xoffsets = 0;

    // m_dataOffsets allows offsets to be precalculated in the constructor then simply looked up according to loop
    // indices during runtime.
    constexpr static std::array<int32, TP_INTERPOLATE_FACTOR> m_dataStarts =
        fnInitStarts<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF>();
    constexpr static std::array<uint64, TP_INTERPOLATE_FACTOR> m_dataOffsets =
        fnInitOffsets<TP_INTERPOLATE_FACTOR, TT_DATA, TT_COEFF>();
    window_incr(inInterface.inWindow, m_kDataWindowOffset); // move input data pointer past the margin padding

// preamble, load data from window into register
// load the data registers with enough data for the initial MAC(MUL)
#pragma unroll(m_kInitialLoadsIncr - 1)
    for (int initLoads = 0; initLoads < m_kInitialLoadsIncr - 1; ++initLoads) {
        readData = window_readincr_256b<TT_DATA>(inInterface.inWindow);
        sbuff.val = upd_w(sbuff.val, initLoads, readData.val);
    }

    // loop through window, computing a vector of output for each iteration.
    for (unsigned i = 0; i < m_kLsize / m_kRepeatFactor; i++)
        chess_prepare_for_pipelining chess_loop_range(m_kLsize / m_kRepeatFactor, ) {
            dataLoaded = 0;
            dataNeeded = m_kVOutSize;

            numDataLoads = m_kInitialLoadsIncr - 1;
#pragma unroll(m_kRepeatFactor)
            for (int strobe = 0; strobe < m_kRepeatFactor; strobe++) {
                if (dataNeeded > dataLoaded) {
                    readData = window_readincr_256b<TT_DATA>(inInterface.inWindow);
                    sbuff.val = upd_w(sbuff.val, numDataLoads % m_kDataLoadsInReg, readData.val);
                    dataLoaded += m_kDataLoadVsize;
                    numDataLoads++;
                }
                dataNeeded += m_kLanes;
#pragma unroll(m_kPhases)
                // The phase loop effectively multiplies the number of lanes in use to ensures that
                // an integer number of interpolation polyphases are calculated
                for (int phase = 0; phase < m_kPhases; ++phase) {
                    coeff = ((T_buff_256b<TT_COEFF>*)&m_internalTaps[LCMPhase][0][0][0]);
                    coe0 = *coeff;
                    xstart = m_kDataBuffXOffset + m_dataStarts[phase] + strobe * m_kLanes;

                    // Read cascade input. Do nothing if cascade input not present.
                    acc = readCascade<TT_DATA, TT_COEFF>(inInterface, acc);
                    // Perform the first term (or 2 or 4) of the FIR polynomial.
                    acc = initMacIntAsym(inInterface, acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes,
                                         m_dataOffsets[LCMPhase], xstart);

// loop to work through the operations in the FIR polynomial
#pragma unroll(GUARD_ZERO((m_kNumOps - 1)))
                    for (int op = 1; op < m_kNumOps; ++op) {
                        coeff = ((T_buff_256b<TT_COEFF>*)&m_internalTaps[LCMPhase][op][0][0]);
                        coe0 = *coeff;
                        xstart += m_kColumns;
                        // perform the MAC, i.e. cacculate some terms of the FIR polynomial
                        acc = macIntAsym(acc, sbuff, coe0, TP_INTERPOLATE_FACTOR, m_kLanes, m_dataOffsets[LCMPhase],
                                         xstart);
                    }

                    // Write cascade. Do nothing if cascade not present.
                    writeCascade<TT_DATA, TT_COEFF>(outInterface, acc);

                    // convert the accumulator type into an integer type for output, downshift and apply any rounding
                    // and saturation.
                    outVal = shiftAndSaturateIntAsym(acc, TP_SHIFT);
                    // Write to output window
                    writeWindow<TT_DATA, TT_COEFF, TP_NUM_OUTPUTS>(outInterface, outVal);

                    LCMPhase = (LCMPhase + 1);
                    if (LCMPhase == m_kLCMPhases) { // a mod function without the division.TODO optimize TP_INTERPOLATE
                                                    // to prime factor
                        LCMPhase = 0;
                    }
                }
            } // strobe loop
        }     // i loop (Lsize)
};

// FIR filter function overloaded with cascade interface variations
// Default declaration and specialization for one kernel, static coefficients, single output
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          bool TP_CASC_IN,
          bool TP_CASC_OUT,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_USE_COEFF_RELOAD,
          unsigned int TP_DUAL_IP,
          unsigned int TP_NUM_OUTPUTS,
          unsigned int TP_API,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          TP_CASC_IN,
                          TP_CASC_OUT,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          TP_USE_COEFF_RELOAD,
                          TP_DUAL_IP,
                          TP_NUM_OUTPUTS,
                          TP_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       output_window<TT_DATA>* outWindow) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    outInterface.outWindow = outWindow;
    this->filterKernel(inInterface, outInterface);
};

// Static coefficients, dual output
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       output_window<TT_DATA>* outWindow,
                                                       output_window<TT_DATA>* outWindow2) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    outInterface.outWindow = outWindow;
    outInterface.outWindow2 = outWindow2;
    this->filterKernel(inInterface, outInterface);
};

// This is a specialization of the main class for when there is only one kernel for the whole filter.
// Reloadable coefficients, single output
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       output_window<TT_DATA>* outWindow,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    outInterface.outWindow = outWindow;
    this->filterKernel(inInterface, outInterface, inTaps);
};

// Reloadable coefficients, single output
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       output_window<TT_DATA>* outWindow,
                                                       output_window<TT_DATA>* outWindow2,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    outInterface.outWindow = outWindow;
    outInterface.outWindow2 = outWindow2;
    this->filterKernel(inInterface, outInterface, inTaps);
};

// This is a specialization of the main class for the final kernel in a cascade chain.
// Static coefficients, single output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       input_stream_cacc48* inCascade,
                                                       output_window<TT_DATA>* outWindow) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    inInterface.inCascade = inCascade;
    outInterface.outWindow = outWindow;
    this->filterKernel(inInterface, outInterface);
};

// Static coefficients, dual output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       input_stream_cacc48* inCascade,
                                                       output_window<TT_DATA>* outWindow,
                                                       output_window<TT_DATA>* outWindow2) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    inInterface.inCascade = inCascade;
    outInterface.outWindow = outWindow;
    outInterface.outWindow2 = outWindow2;
    this->filterKernel(inInterface, outInterface);
};

// FIR filter function overloaded with cascade interface variations
// This is a specialization of the main class for the first kernel in a cascade chain.
// Static coefficients
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       output_stream_cacc48* outCascade,
                                                       output_window<TT_DATA>* broadcastWindow) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    outInterface.outCascade = outCascade;
    outInterface.broadcastWindow = broadcastWindow;
    this->filterKernel(inInterface, outInterface);
};

// FIR filter function overloaded with cascade interface variations
// This is a specialization of the main class for any kernel within a cascade chain, but neither first nor last.
// Static coefficients
//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream_cacc48* outCascade,
                                                       output_window<TT_DATA>* broadcastWindow) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    inInterface.inCascade = inCascade;
    outInterface.outCascade = outCascade;
    outInterface.broadcastWindow = broadcastWindow;
    this->filterKernel(inInterface, outInterface);
};

// This is a specialization of the main class for the final kernel in a cascade chain.
// Reloadable coefficients, single output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       input_stream_cacc48* inCascade,
                                                       output_window<TT_DATA>* outWindow) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    inInterface.inCascade = inCascade;
    outInterface.outWindow = outWindow;
    this->filterKernelRtp(inInterface, outInterface);
};

// This is a specialization of the main class for the final kernel in a cascade chain.
// Reloadable coefficients, dual output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       input_stream_cacc48* inCascade,
                                                       output_window<TT_DATA>* outWindow,
                                                       output_window<TT_DATA>* outWindow2) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    inInterface.inCascade = inCascade;
    outInterface.outWindow = outWindow;
    outInterface.outWindow2 = outWindow2;
    this->filterKernelRtp(inInterface, outInterface);
};

// This is a specialization of the main class for the first kernel in a cascade chain.
// Reloadable coefficients
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       output_stream_cacc48* outCascade,
                                                       output_window<TT_DATA>* broadcastWindow,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    outInterface.outCascade = outCascade;
    outInterface.broadcastWindow = broadcastWindow;
    this->filterKernel(inInterface, outInterface, inTaps);
};

// This is a specialization of the main class for any kernel within a cascade chain, but neither first nor last.
// Reloadable coefficients
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_WINDOW_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_window<TT_DATA>* inWindow,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream_cacc48* outCascade,
                                                       output_window<TT_DATA>* broadcastWindow) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inWindow = inWindow;
    inInterface.inCascade = inCascade;
    outInterface.outCascade = outCascade;
    outInterface.broadcastWindow = broadcastWindow;
    this->filterKernelRtp(inInterface, outInterface);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////  STREAM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  Single Kernel. Static Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       output_stream<TT_DATA>* outStream) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    outInterface.outStream = outStream;
    this->filterKernel(inInterface, outInterface);
};

//  Single Kernel. Static Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernel(inInterface, outInterface);
};

//  Single Kernel. Reloadable Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       output_stream<TT_DATA>* outStream,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    outInterface.outStream = outStream;
    this->filterKernel(inInterface, outInterface, inTaps);
};

//  Single Kernel. Reloadable Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernel(inInterface, outInterface, inTaps);
};

//  Cascaded Kernels - Final Kernel. Static Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - Final Kernel. Static Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - First Kernel. Static Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       output_stream_cacc48* outCascade) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    outInterface.outCascade = outCascade;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - Middle Kernel. Static Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream_cacc48* outCascade) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inCascade = inCascade;
    outInterface.outCascade = outCascade;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - Final Kernel. Reloadable Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    this->filterKernelRtp(inInterface, outInterface);
};

//  Cascaded Kernels - Final Kernel. Reloadable Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernelRtp(inInterface, outInterface);
};

//  Cascaded Kernels - First Kernel. Reloadable Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       output_stream_cacc48* outCascade,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    outInterface.outCascade = outCascade;
    this->filterKernel(inInterface, outInterface, inTaps);
};

//  Cascaded Kernels - Middle Kernel. Reloadable Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_SINGLE,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream_cacc48* outCascade) {
    T_inputIF<CASC_IN_TRUE, TT_DATA> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inCascade = inCascade;
    outInterface.outCascade = outCascade;
    this->filterKernelRtp(inInterface, outInterface);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Dual Stream

//  Single Kernel. Static Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_DUAL,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       output_stream<TT_DATA>* outStream) {
    T_inputIF<CASC_IN_FALSE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    outInterface.outStream = outStream;
    this->filterKernel(inInterface, outInterface);
};

//  Single Kernel. Static Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_DUAL,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2) {
    T_inputIF<CASC_IN_FALSE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernel(inInterface, outInterface);
};

//  Single Kernel. Reloadable Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_DUAL,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       output_stream<TT_DATA>* outStream,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    outInterface.outStream = outStream;
    this->filterKernel(inInterface, outInterface, inTaps);
};

//  Single Kernel. Reloadable Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          0,
                          1,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_DUAL,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernel(inInterface, outInterface, inTaps);
};

//  Cascaded Kernels - Final Kernel. Static Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_DUAL,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream) {
    T_inputIF<CASC_IN_TRUE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - Final Kernel. Static Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_DUAL,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2) {
    T_inputIF<CASC_IN_TRUE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - First Kernel. Static Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_DUAL,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       output_stream_cacc48* outCascade) {
    T_inputIF<CASC_IN_FALSE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    outInterface.outCascade = outCascade;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - Middle Kernel. Static Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_FALSE,
                          DUAL_IP_DUAL,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream_cacc48* outCascade) {
    T_inputIF<CASC_IN_TRUE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    inInterface.inCascade = inCascade;
    outInterface.outCascade = outCascade;
    this->filterKernel(inInterface, outInterface);
};

//  Cascaded Kernels - Final Kernel. Reloadable Coefficients. Single Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_DUAL,
                          1,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream) {
    T_inputIF<CASC_IN_TRUE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    this->filterKernelRtp(inInterface, outInterface);
};

//  Cascaded Kernels - Final Kernel. Reloadable Coefficients. Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_FALSE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_DUAL,
                          2,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream<TT_DATA>* outStream,
                                                       output_stream<TT_DATA>* outStream2) {
    T_inputIF<CASC_IN_TRUE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_FALSE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    inInterface.inCascade = inCascade;
    outInterface.outStream = outStream;
    outInterface.outStream2 = outStream2;
    this->filterKernelRtp(inInterface, outInterface);
};

//  Cascaded Kernels - First Kernel. Reloadable Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_FALSE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_DUAL,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       output_stream_cacc48* outCascade,
                                                       const TT_COEFF (&inTaps)[TP_COEFF_PHASES_LEN]) {
    T_inputIF<CASC_IN_FALSE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    outInterface.outCascade = outCascade;
    this->filterKernel(inInterface, outInterface, inTaps);
};

//  Cascaded Kernels - Middle Kernel. Reloadable Coefficients. Single/Dual Output
template <typename TT_DATA,
          typename TT_COEFF,
          unsigned int TP_FIR_LEN,
          unsigned int TP_INTERPOLATE_FACTOR,
          unsigned int TP_SHIFT,
          unsigned int TP_RND,
          unsigned int TP_INPUT_WINDOW_VSIZE,
          unsigned int TP_FIR_RANGE_LEN,
          unsigned int TP_KERNEL_POSITION,
          unsigned int TP_CASC_LEN,
          unsigned int TP_NUM_OUTPUTS,
          int TP_MODIFY_MARGIN_OFFSET,
          unsigned int TP_COEFF_PHASE,
          unsigned int TP_COEFF_PHASE_OFFSET,
          unsigned int TP_COEFF_PHASES,
          unsigned int TP_COEFF_PHASES_LEN>
void fir_interpolate_asym<TT_DATA,
                          TT_COEFF,
                          TP_FIR_LEN,
                          TP_INTERPOLATE_FACTOR,
                          TP_SHIFT,
                          TP_RND,
                          TP_INPUT_WINDOW_VSIZE,
                          CASC_IN_TRUE,
                          CASC_OUT_TRUE,
                          TP_FIR_RANGE_LEN,
                          TP_KERNEL_POSITION,
                          TP_CASC_LEN,
                          USE_COEFF_RELOAD_TRUE,
                          DUAL_IP_DUAL,
                          TP_NUM_OUTPUTS,
                          USE_STREAM_API,
                          TP_MODIFY_MARGIN_OFFSET,
                          TP_COEFF_PHASE,
                          TP_COEFF_PHASE_OFFSET,
                          TP_COEFF_PHASES,
                          TP_COEFF_PHASES_LEN>::filter(input_stream<TT_DATA>* inStream,
                                                       input_stream<TT_DATA>* inStream2,
                                                       input_stream_cacc48* inCascade,
                                                       output_stream_cacc48* outCascade) {
    T_inputIF<CASC_IN_TRUE, TT_DATA, DUAL_IP_DUAL> inInterface;
    T_outputIF<CASC_OUT_TRUE, TT_DATA> outInterface;
    inInterface.inStream = inStream;
    inInterface.inStream2 = inStream2;
    inInterface.inCascade = inCascade;
    outInterface.outCascade = outCascade;
    this->filterKernelRtp(inInterface, outInterface);
};
}
}
}
}
} // namespaces
