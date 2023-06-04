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
//------------------------------------------------------------------------------
// UUT DEFAULT CONFIGURATION
#ifndef DATA_TYPE
#define DATA_TYPE cint16
#endif
#ifndef COEFF_TYPE
#define COEFF_TYPE int16
#endif
#ifndef FIR_LEN
#define FIR_LEN 32
#endif
#ifndef INTERPOLATE_FACTOR
#define INTERPOLATE_FACTOR 3
#endif
#ifndef DECIMATE_FACTOR
#define DECIMATE_FACTOR 2
#endif
#ifndef SHIFT
#define SHIFT 16
#endif
#ifndef ROUND_MODE
#define ROUND_MODE 0
#endif
#ifndef INPUT_WINDOW_VSIZE
#define INPUT_WINDOW_VSIZE 256
#endif
#ifndef USE_CHAIN
#define USE_CHAIN 0
#endif

#ifndef INPUT_FILE
#define INPUT_FILE "data/input.txt"
#endif
#ifndef OUTPUT_FILE
#define OUTPUT_FILE "data/output.txt"
#endif

#ifndef NUM_ITER
#define NUM_ITER 1
#endif

#ifndef PORT_API
#define PORT_API 0
#endif

#ifndef UUT_SSR
#define UUT_SSR 1
#endif

#ifndef CASC_LEN
#define CASC_LEN 1
#endif

#ifndef USE_COEFF_RELOAD
#define USE_COEFF_RELOAD 0
#endif

#define INPUT_SAMPLES INPUT_WINDOW_VSIZE* NUM_ITER
#define INPUT_MARGIN(x, y) CEIL(x, (32 / sizeof(y)))
#define OUTPUT_SAMPLES (INPUT_WINDOW_VSIZE * NUM_ITER * INTERPOLATE_FACTOR) / DECIMATE_FACTOR

#ifndef COEFF_SEED
#define COEFF_SEED 0xC0FFEE
#endif

// Force reference model to ignore SSR and only use array length 1
#ifdef USING_UUT
#define P_SSR UUT_SSR
#else
#define P_SSR 1
#endif

#ifndef USING_UUT
#undef PORT_API
#undef NUM_OUTPUTS
#undef DUAL_IP
#define PORT_API 0
#define NUM_OUTPUTS 1
#define DUAL_IP 0
#endif

#ifndef DUAL_IP
#define DUAL_IP 0
#endif

// END OF UUT CONFIGURATION
//------------------------------------------------------------------------------
