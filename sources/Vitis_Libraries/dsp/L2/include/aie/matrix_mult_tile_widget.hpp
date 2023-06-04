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
#ifndef _DSPLIB_CONDITIONAL_WIDGET_HPP_
#define _DSPLIB_CONDITIONAL_WIDGET_HPP_

// This file holds the definition of the conditional widget class
/**
 * @file conditional_widget.hpp
 *
 **/

#include <adf.h>
#include <vector>

namespace xf {
namespace dsp {
namespace aie {
namespace blas {
namespace matrix_mult {

using namespace adf;

/**
 * @brief ConditionalWidget conditionally instances a kernel of a defined type.
 * In addition, ConditionalWidget will connect graph's input and output ports to kernel's ports. Ports must be of window
 *type.
 * When disabled, no kernel will be created and input wil be directly connected to output.
 *
 * These are the templates to configure the ConditionalWidget graph class.
 *
 * @ingroup gemm_graph
 *
 * @tparam addWidget conditionally add a widget kernel.  \n
 *         When set to 1, kernel of: ``class widgetClass`` will be created.
 * @tparam windowSize describes the size of the window the ``widgetClass`` is to operate on.
 * @tparam widgetClass defines the class to be created.
 *
**/
template <unsigned int addWidget, unsigned int windowSize, class widgetClass>
class ConditionalWidget {
   public:
    using portConnect = connect<window<windowSize> >;
    ConditionalWidget(){}; // default constructor
    template <typename inType, typename outType>
    static kernel create(port<inType>& inPort, port<outType>& outPort) {
        kernel widget;
        if (addWidget == 1) {
            widget = kernel::create_object<widgetClass>();
            portConnect(inPort, widget.in[0]);
            portConnect(widget.out[0], outPort);
        } else {
            portConnect(inPort, outPort);
        }

        return widget;
    }
};
}
}
}
}
}

#endif