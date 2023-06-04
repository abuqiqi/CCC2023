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
#ifndef _DSPLIB_widget_api_cast_REF_HPP_
#define _DSPLIB_widget_api_cast_REF_HPP_

/*
Widget API Cast reference model
*/

#include <adf.h>
#include <limits>
#include "fir_ref_utils.hpp"
#include "widget_api_cast_traits.hpp"

#ifndef _DSPLIB_WIDGET_API_CAST_REF_DEBUG_
//#define _DSPLIB_WIDGET_API_CAST_REF_DEBUG_
//#include "debug_utils.h"
#endif //_DSPLIB_WIDGET_API_CAST_REF_DEBUG_

namespace xf {
namespace dsp {
namespace aie {
namespace widget {
namespace api_cast {

//-----------------------------------------------------------------------------------------------------
template <typename TT_DATA, // type of data input and output
          unsigned int TP_IN_API,
          unsigned int TP_OUT_API,
          unsigned int TP_NUM_INPUTS,
          unsigned int TP_WINDOW_VSIZE, // does not include header, so refers to payload samples
          unsigned int TP_NUM_OUTPUT_CLONES = 1,
          unsigned int TP_PATTERN = 0,
          unsigned int TP_HEADER_BYTES = 0>
class widget_api_cast_ref {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_window<TT_DATA>* inWindow0, output_window<TT_DATA>* outWindow0);
};

// window to window, 1 to 2
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kWindowAPI, kWindowAPI, 1, TP_WINDOW_VSIZE, 2, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_window<TT_DATA>* inWindow0,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1);
};

// window to window, 1 to 3
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kWindowAPI, kWindowAPI, 1, TP_WINDOW_VSIZE, 3, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_window<TT_DATA>* inWindow0,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1,
                      output_window<TT_DATA>* outWindow2);
};

// stream to window, 1 to 1
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 1, TP_WINDOW_VSIZE, 1, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0, output_window<TT_DATA>* outWindow0);
};

// stream to window, 1 to 2
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 1, TP_WINDOW_VSIZE, 2, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1);
};

// stream to window, 1 to 3
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 1, TP_WINDOW_VSIZE, 3, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1,
                      output_window<TT_DATA>* outWindow2);
};

// stream to window, 1 to 4
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 1, TP_WINDOW_VSIZE, 4, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }

    void transferData(input_stream<TT_DATA>* inStream0,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1,
                      output_window<TT_DATA>* outWindow2,
                      output_window<TT_DATA>* outWindow3);
};

// Dual stream in
// stream to window, 2 to 1
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 2, TP_WINDOW_VSIZE, 1, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0,
                      input_stream<TT_DATA>* inStream1,
                      output_window<TT_DATA>* outWindow0);
};

// stream to window, 2 to 2
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 2, TP_WINDOW_VSIZE, 2, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0,
                      input_stream<TT_DATA>* inStream1,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1);
};

// stream to window, 2 to 3
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 2, TP_WINDOW_VSIZE, 3, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0,
                      input_stream<TT_DATA>* inStream1,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1,
                      output_window<TT_DATA>* outWindow2);
};

// stream to window, 2 to 4
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kStreamAPI, kWindowAPI, 2, TP_WINDOW_VSIZE, 4, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_stream<TT_DATA>* inStream0,
                      input_stream<TT_DATA>* inStream1,
                      output_window<TT_DATA>* outWindow0,
                      output_window<TT_DATA>* outWindow1,
                      output_window<TT_DATA>* outWindow2,
                      output_window<TT_DATA>* outWindow3);
};

// Window to stream
// window to stream, 1 to 1
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kWindowAPI, kStreamAPI, 1, TP_WINDOW_VSIZE, 1, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_window<TT_DATA>* inWindow0, output_stream<TT_DATA>* outStream0);
};

// window to stream, 1 to 2
template <typename TT_DATA, // type of data input and output
          unsigned int TP_WINDOW_VSIZE,
          unsigned int TP_PATTERN,
          unsigned int TP_HEADER_BYTES>
class widget_api_cast_ref<TT_DATA, kWindowAPI, kStreamAPI, 1, TP_WINDOW_VSIZE, 2, TP_PATTERN, TP_HEADER_BYTES> {
   private:
   public:
    int kIndex;
    widget_api_cast_ref(int idx) { kIndex = idx; }

    // Register Kernel Class
    static void registerKernelClass() { REGISTER_FUNCTION(widget_api_cast_ref::transferData); }
    void transferData(input_window<TT_DATA>* inWindow0,
                      output_stream<TT_DATA>* outStream0,
                      output_stream<TT_DATA>* outStream1);
};
}
}
}
}
}

#endif // _DSPLIB_widget_api_cast_REF_HPP_
