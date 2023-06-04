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
This file holds the body of the test harness for dds_mixer reference model graph
*/

#include <stdio.h>
#include "test.hpp"

xf::dsp::aie::testcase::test_graph ddsMix;

int main(void) {
    printf("\n");
    printf("========================\n");
    printf("TEST.CPP STARTED\n");
    printf("UUT: ");
    printf(QUOTE(UUT_GRAPH));
    printf("\n");
    printf("========================\n");
    printf("Input stimulus file");
    printf(QUOTE(INPUT_FILE));
    printf("\n");
    printf("Input stimulus file 2");
    printf(QUOTE(INPUT_FILE2));
    printf("\n");
    printf("Input samples   = %d \n", INPUT_SAMPLES);
    printf("Output samples  = %d \n", OUTPUT_SAMPLES);
    printf("Data type       = ");
    printf(QUOTE(DATA_TYPE));
    printf("\n");
    printf("TEST.CPP PRINT FINISHED \n");
    printf("\n");

    ddsMix.init();

    ddsMix.run(NITER);

    ddsMix.end();

    printf("TEST.CPP IS FINISHED \n");
    printf("========================\n");
    printf("========================\n");

    return 0;
}
