#pragma once

#include <aie_api/aie.hpp>
#include <aie_api/aie_adf.hpp>

#define N_POINT 1024
#define MAX_VEC_LEN 8
#define N_SWAP_2 240
#define N_SWAP_4 768
#define MAT_OMG_SHIFT 14
#define OMG_SHIFT 14

using namespace aie;

const int16 swap2[]={2, 256, 5, 640, 7, 896, 8, 64, 10, 320, 13, 704, 15, 960, 16, 32, 18, 288, 21, 672, 23, 928, 24, 96, 26, 352, 29, 736, 31, 992, 34, 272, 37, 656, 39, 912, 40, 80, 42, 336, 45, 720, 47, 976, 50, 304, 53, 688, 55, 944, 56, 112, 58, 368, 61, 752, 63, 1008, 66, 264, 69, 648, 71, 904, 74, 328, 77, 712, 79, 968, 82, 296, 85, 680, 87, 936, 88, 104, 90, 360, 93, 744, 95, 1000, 98, 280, 101, 664, 103, 920, 106, 344, 109, 728, 111, 984, 114, 312, 117, 696, 119, 952, 122, 376, 125, 760, 127, 1016, 261, 642, 263, 898, 266, 322, 269, 706, 271, 962, 274, 290, 277, 674, 279, 930, 282, 354, 285, 738, 287, 994, 293, 658, 295, 914, 298, 338, 301, 722, 303, 978, 309, 690, 311, 946, 314, 370, 317, 754, 319, 1010, 325, 650, 327, 906, 333, 714, 335, 970, 341, 682, 343, 938, 346, 362, 349, 746, 351, 1002, 357, 666, 359, 922, 365, 730, 367, 986, 373, 698, 375, 954, 381, 762, 383, 1018, 647, 901, 653, 709, 655, 965, 661, 677, 663, 933, 669, 741, 671, 997, 679, 917, 685, 725, 687, 981, 695, 949, 701, 757, 703, 1013, 711, 909, 719, 973, 727, 941, 733, 749, 735, 1005, 743, 925, 751, 989, 759, 957, 767, 1021, 911, 967, 919, 935, 927, 999, 943, 983, 959, 1015, 991, 1007};
const int16 swap4[]={1, 128, 4, 512, 3, 384, 6, 768, 9, 192, 12, 576, 11, 448, 14, 832, 17, 160, 20, 544, 19, 416, 22, 800, 25, 224, 28, 608, 27, 480, 30, 864, 33, 144, 36, 528, 35, 400, 38, 784, 41, 208, 44, 592, 43, 464, 46, 848, 49, 176, 52, 560, 51, 432, 54, 816, 57, 240, 60, 624, 59, 496, 62, 880, 65, 136, 68, 520, 67, 392, 70, 776, 73, 200, 76, 584, 75, 456, 78, 840, 81, 168, 84, 552, 83, 424, 86, 808, 89, 232, 92, 616, 91, 488, 94, 872, 97, 152, 100, 536, 99, 408, 102, 792, 105, 216, 108, 600, 107, 472, 110, 856, 113, 184, 116, 568, 115, 440, 118, 824, 121, 248, 124, 632, 123, 504, 126, 888, 129, 132, 516, 513, 130, 260, 514, 257, 131, 388, 518, 769, 133, 644, 517, 641, 134, 772, 515, 385, 135, 900, 519, 897, 137, 196, 524, 577, 138, 324, 522, 321, 139, 452, 526, 833, 140, 580, 521, 193, 141, 708, 525, 705, 142, 836, 523, 449, 143, 964, 527, 961, 145, 164, 532, 545, 146, 292, 530, 289, 147, 420, 534, 801, 148, 548, 529, 161, 149, 676, 533, 673, 150, 804, 531, 417, 151, 932, 535, 929, 153, 228, 540, 609, 154, 356, 538, 353, 155, 484, 542, 865, 156, 612, 537, 225, 157, 740, 541, 737, 158, 868, 539, 481, 159, 996, 543, 993, 162, 276, 546, 273, 163, 404, 550, 785, 165, 660, 549, 657, 166, 788, 547, 401, 167, 916, 551, 913, 169, 212, 556, 593, 170, 340, 554, 337, 171, 468, 558, 849, 172, 596, 553, 209, 173, 724, 557, 721, 174, 852, 555, 465, 175, 980, 559, 977, 177, 180, 564, 561, 178, 308, 562, 305, 179, 436, 566, 817, 181, 692, 565, 689, 182, 820, 563, 433, 183, 948, 567, 945, 185, 244, 572, 625, 186, 372, 570, 369, 187, 500, 574, 881, 188, 628, 569, 241, 189, 756, 573, 753, 190, 884, 571, 497, 191, 1012, 575, 1009, 194, 268, 578, 265, 195, 396, 582, 777, 197, 652, 581, 649, 198, 780, 579, 393, 199, 908, 583, 905, 201, 204, 588, 585, 202, 332, 586, 329, 203, 460, 590, 841, 205, 716, 589, 713, 206, 844, 587, 457, 207, 972, 591, 969, 210, 300, 594, 297, 211, 428, 598, 809, 213, 684, 597, 681, 214, 812, 595, 425, 215, 940, 599, 937, 217, 236, 604, 617, 218, 364, 602, 361, 219, 492, 606, 873, 220, 620, 601, 233, 221, 748, 605, 745, 222, 876, 603, 489, 223, 1004, 607, 1001, 226, 284, 610, 281, 227, 412, 614, 793, 229, 668, 613, 665, 230, 796, 611, 409, 231, 924, 615, 921, 234, 348, 618, 345, 235, 476, 622, 857, 237, 732, 621, 729, 238, 860, 619, 473, 239, 988, 623, 985, 242, 316, 626, 313, 243, 444, 630, 825, 245, 700, 629, 697, 246, 828, 627, 441, 247, 956, 631, 953, 249, 252, 636, 633, 250, 380, 634, 377, 251, 508, 638, 889, 253, 764, 637, 761, 254, 892, 635, 505, 255, 1020, 639, 1017, 259, 386, 262, 770, 267, 450, 270, 834, 275, 418, 278, 802, 283, 482, 286, 866, 291, 402, 294, 786, 299, 466, 302, 850, 307, 434, 310, 818, 315, 498, 318, 882, 323, 394, 326, 778, 331, 458, 334, 842, 339, 426, 342, 810, 347, 490, 350, 874, 355, 410, 358, 794, 363, 474, 366, 858, 371, 442, 374, 826, 379, 506, 382, 890, 387, 390, 774, 771, 389, 646, 773, 643, 391, 902, 775, 899, 395, 454, 782, 835, 397, 710, 781, 707, 398, 838, 779, 451, 399, 966, 783, 963, 403, 422, 790, 803, 405, 678, 789, 675, 406, 806, 787, 419, 407, 934, 791, 931, 411, 486, 798, 867, 413, 742, 797, 739, 414, 870, 795, 483, 415, 998, 799, 995, 421, 662, 805, 659, 423, 918, 807, 915, 427, 470, 814, 851, 429, 726, 813, 723, 430, 854, 811, 467, 431, 982, 815, 979, 435, 438, 822, 819, 437, 694, 821, 691, 439, 950, 823, 947, 443, 502, 830, 883, 445, 758, 829, 755, 446, 886, 827, 499, 447, 1014, 831, 1011, 453, 654, 837, 651, 455, 910, 839, 907, 459, 462, 846, 843, 461, 718, 845, 715, 463, 974, 847, 971, 469, 686, 853, 683, 471, 942, 855, 939, 475, 494, 862, 875, 477, 750, 861, 747, 478, 878, 859, 491, 479, 1006, 863, 1003, 485, 670, 869, 667, 487, 926, 871, 923, 493, 734, 877, 731, 495, 990, 879, 987, 501, 702, 885, 699, 503, 958, 887, 955, 507, 510, 894, 891, 509, 766, 893, 763, 511, 1022, 895, 1019};

static cint16 mat_omg_8[64]={{16384,0},{16384,0},{16384,0},{16384,0},{16384,0},{16384,0},{16384,0},{16384,0},{16384,0},{11585,-11585},{0,-16384},{-11585,-11585},{-16384,0},{-11585,11585},{0,16384},{11585,11585},{16384,0},{0,-16384},{-16384,0},{0,16384},{16384,0},{0,-16384},{-16384,0},{0,16384},{16384,0},{-11585,-11585},{0,16384},{11585,-11585},{-16384,0},{11585,11585},{0,-16384},{-11585,11585},{16384,0},{-16384,0},{16384,0},{-16384,0},{16384,0},{-16384,0},{16384,0},{-16384,0},{16384,0},{-11585,11585},{0,-16384},{11585,11585},{-16384,0},{11585,-11585},{0,16384},{-11585,-11585},{16384,0},{0,16384},{-16384,0},{0,-16384},{16384,0},{0,16384},{-16384,0},{0,-16384},{16384,0},{11585,11585},{0,16384},{-11585,11585},{-16384,0},{-11585,-11585},{0,-16384},{11585,-11585},};

static cint16 omg_16[]={{16384,0},{15136,-6269},{11585,-11585},{6269,-15136},{0,-16384},{-6269,-15136},{-11585,-11585},{-15136,-6269},};

static cint16 omg_32[]={{16384,0},{16069,-3196},{15136,-6269},{13622,-9102},{11585,-11585},{9102,-13622},{6269,-15136},{3196,-16069},{0,-16384},{-3196,-16069},{-6269,-15136},{-9102,-13622},{-11585,-11585},{-13622,-9102},{-15136,-6269},{-16069,-3196},};

static cint16 omg_64[]={{16384,0},{16305,-1605},{16069,-3196},{15678,-4756},{15136,-6269},{14449,-7723},{13622,-9102},{12665,-10393},{11585,-11585},{10393,-12665},{9102,-13622},{7723,-14449},{6269,-15136},{4756,-15678},{3196,-16069},{1605,-16305},{0,-16384},{-1605,-16305},{-3196,-16069},{-4756,-15678},{-6269,-15136},{-7723,-14449},{-9102,-13622},{-10393,-12665},{-11585,-11585},{-12665,-10393},{-13622,-9102},{-14449,-7723},{-15136,-6269},{-15678,-4756},{-16069,-3196},{-16305,-1605},};

static cint16 omg_128[]={{16384,0},{16364,-803},{16305,-1605},{16206,-2404},{16069,-3196},{15892,-3980},{15678,-4756},{15426,-5519},{15136,-6269},{14810,-7005},{14449,-7723},{14053,-8423},{13622,-9102},{13159,-9759},{12665,-10393},{12139,-11002},{11585,-11585},{11002,-12139},{10393,-12665},{9759,-13159},{9102,-13622},{8423,-14053},{7723,-14449},{7005,-14810},{6269,-15136},{5519,-15426},{4756,-15678},{3980,-15892},{3196,-16069},{2404,-16206},{1605,-16305},{803,-16364},{0,-16384},{-803,-16364},{-1605,-16305},{-2404,-16206},{-3196,-16069},{-3980,-15892},{-4756,-15678},{-5519,-15426},{-6269,-15136},{-7005,-14810},{-7723,-14449},{-8423,-14053},{-9102,-13622},{-9759,-13159},{-10393,-12665},{-11002,-12139},{-11585,-11585},{-12139,-11002},{-12665,-10393},{-13159,-9759},{-13622,-9102},{-14053,-8423},{-14449,-7723},{-14810,-7005},{-15136,-6269},{-15426,-5519},{-15678,-4756},{-15892,-3980},{-16069,-3196},{-16206,-2404},{-16305,-1605},{-16364,-803},};

static cint16 omg_256[]={{16384,0},{16379,-402},{16364,-803},{16339,-1205},{16305,-1605},{16260,-2005},{16206,-2404},{16142,-2801},{16069,-3196},{15985,-3589},{15892,-3980},{15790,-4369},{15678,-4756},{15557,-5139},{15426,-5519},{15286,-5896},{15136,-6269},{14978,-6639},{14810,-7005},{14634,-7366},{14449,-7723},{14255,-8075},{14053,-8423},{13842,-8765},{13622,-9102},{13395,-9434},{13159,-9759},{12916,-10079},{12665,-10393},{12406,-10701},{12139,-11002},{11866,-11297},{11585,-11585},{11297,-11866},{11002,-12139},{10701,-12406},{10393,-12665},{10079,-12916},{9759,-13159},{9434,-13395},{9102,-13622},{8765,-13842},{8423,-14053},{8075,-14255},{7723,-14449},{7366,-14634},{7005,-14810},{6639,-14978},{6269,-15136},{5896,-15286},{5519,-15426},{5139,-15557},{4756,-15678},{4369,-15790},{3980,-15892},{3589,-15985},{3196,-16069},{2801,-16142},{2404,-16206},{2005,-16260},{1605,-16305},{1205,-16339},{803,-16364},{402,-16379},{0,-16384},{-402,-16379},{-803,-16364},{-1205,-16339},{-1605,-16305},{-2005,-16260},{-2404,-16206},{-2801,-16142},{-3196,-16069},{-3589,-15985},{-3980,-15892},{-4369,-15790},{-4756,-15678},{-5139,-15557},{-5519,-15426},{-5896,-15286},{-6269,-15136},{-6639,-14978},{-7005,-14810},{-7366,-14634},{-7723,-14449},{-8075,-14255},{-8423,-14053},{-8765,-13842},{-9102,-13622},{-9434,-13395},{-9759,-13159},{-10079,-12916},{-10393,-12665},{-10701,-12406},{-11002,-12139},{-11297,-11866},{-11585,-11585},{-11866,-11297},{-12139,-11002},{-12406,-10701},{-12665,-10393},{-12916,-10079},{-13159,-9759},{-13395,-9434},{-13622,-9102},{-13842,-8765},{-14053,-8423},{-14255,-8075},{-14449,-7723},{-14634,-7366},{-14810,-7005},{-14978,-6639},{-15136,-6269},{-15286,-5896},{-15426,-5519},{-15557,-5139},{-15678,-4756},{-15790,-4369},{-15892,-3980},{-15985,-3589},{-16069,-3196},{-16142,-2801},{-16206,-2404},{-16260,-2005},{-16305,-1605},{-16339,-1205},{-16364,-803},{-16379,-402},};

static cint16 omg_512[]={{16384,0},{16382,-201},{16379,-402},{16372,-603},{16364,-803},{16353,-1004},{16339,-1205},{16323,-1405},{16305,-1605},{16284,-1805},{16260,-2005},{16234,-2204},{16206,-2404},{16175,-2602},{16142,-2801},{16107,-2998},{16069,-3196},{16028,-3393},{15985,-3589},{15940,-3785},{15892,-3980},{15842,-4175},{15790,-4369},{15735,-4563},{15678,-4756},{15618,-4948},{15557,-5139},{15492,-5329},{15426,-5519},{15357,-5708},{15286,-5896},{15212,-6083},{15136,-6269},{15058,-6455},{14978,-6639},{14895,-6822},{14810,-7005},{14723,-7186},{14634,-7366},{14543,-7545},{14449,-7723},{14353,-7900},{14255,-8075},{14155,-8249},{14053,-8423},{13948,-8594},{13842,-8765},{13733,-8934},{13622,-9102},{13510,-9268},{13395,-9434},{13278,-9597},{13159,-9759},{13038,-9920},{12916,-10079},{12791,-10237},{12665,-10393},{12536,-10548},{12406,-10701},{12273,-10853},{12139,-11002},{12003,-11150},{11866,-11297},{11726,-11442},{11585,-11585},{11442,-11726},{11297,-11866},{11150,-12003},{11002,-12139},{10853,-12273},{10701,-12406},{10548,-12536},{10393,-12665},{10237,-12791},{10079,-12916},{9920,-13038},{9759,-13159},{9597,-13278},{9434,-13395},{9268,-13510},{9102,-13622},{8934,-13733},{8765,-13842},{8594,-13948},{8423,-14053},{8249,-14155},{8075,-14255},{7900,-14353},{7723,-14449},{7545,-14543},{7366,-14634},{7186,-14723},{7005,-14810},{6822,-14895},{6639,-14978},{6455,-15058},{6269,-15136},{6083,-15212},{5896,-15286},{5708,-15357},{5519,-15426},{5329,-15492},{5139,-15557},{4948,-15618},{4756,-15678},{4563,-15735},{4369,-15790},{4175,-15842},{3980,-15892},{3785,-15940},{3589,-15985},{3393,-16028},{3196,-16069},{2998,-16107},{2801,-16142},{2602,-16175},{2404,-16206},{2204,-16234},{2005,-16260},{1805,-16284},{1605,-16305},{1405,-16323},{1205,-16339},{1004,-16353},{803,-16364},{603,-16372},{402,-16379},{201,-16382},{0,-16384},{-201,-16382},{-402,-16379},{-603,-16372},{-803,-16364},{-1004,-16353},{-1205,-16339},{-1405,-16323},{-1605,-16305},{-1805,-16284},{-2005,-16260},{-2204,-16234},{-2404,-16206},{-2602,-16175},{-2801,-16142},{-2998,-16107},{-3196,-16069},{-3393,-16028},{-3589,-15985},{-3785,-15940},{-3980,-15892},{-4175,-15842},{-4369,-15790},{-4563,-15735},{-4756,-15678},{-4948,-15618},{-5139,-15557},{-5329,-15492},{-5519,-15426},{-5708,-15357},{-5896,-15286},{-6083,-15212},{-6269,-15136},{-6455,-15058},{-6639,-14978},{-6822,-14895},{-7005,-14810},{-7186,-14723},{-7366,-14634},{-7545,-14543},{-7723,-14449},{-7900,-14353},{-8075,-14255},{-8249,-14155},{-8423,-14053},{-8594,-13948},{-8765,-13842},{-8934,-13733},{-9102,-13622},{-9268,-13510},{-9434,-13395},{-9597,-13278},{-9759,-13159},{-9920,-13038},{-10079,-12916},{-10237,-12791},{-10393,-12665},{-10548,-12536},{-10701,-12406},{-10853,-12273},{-11002,-12139},{-11150,-12003},{-11297,-11866},{-11442,-11726},{-11585,-11585},{-11726,-11442},{-11866,-11297},{-12003,-11150},{-12139,-11002},{-12273,-10853},{-12406,-10701},{-12536,-10548},{-12665,-10393},{-12791,-10237},{-12916,-10079},{-13038,-9920},{-13159,-9759},{-13278,-9597},{-13395,-9434},{-13510,-9268},{-13622,-9102},{-13733,-8934},{-13842,-8765},{-13948,-8594},{-14053,-8423},{-14155,-8249},{-14255,-8075},{-14353,-7900},{-14449,-7723},{-14543,-7545},{-14634,-7366},{-14723,-7186},{-14810,-7005},{-14895,-6822},{-14978,-6639},{-15058,-6455},{-15136,-6269},{-15212,-6083},{-15286,-5896},{-15357,-5708},{-15426,-5519},{-15492,-5329},{-15557,-5139},{-15618,-4948},{-15678,-4756},{-15735,-4563},{-15790,-4369},{-15842,-4175},{-15892,-3980},{-15940,-3785},{-15985,-3589},{-16028,-3393},{-16069,-3196},{-16107,-2998},{-16142,-2801},{-16175,-2602},{-16206,-2404},{-16234,-2204},{-16260,-2005},{-16284,-1805},{-16305,-1605},{-16323,-1405},{-16339,-1205},{-16353,-1004},{-16364,-803},{-16372,-603},{-16379,-402},{-16382,-201},};

static cint16 omg_1024[]={{16384,0},{16383,-100},{16382,-201},{16381,-301},{16379,-402},{16376,-502},{16372,-603},{16368,-703},{16364,-803},{16359,-904},{16353,-1004},{16346,-1105},{16339,-1205},{16331,-1305},{16323,-1405},{16314,-1505},{16305,-1605},{16294,-1705},{16284,-1805},{16272,-1905},{16260,-2005},{16248,-2105},{16234,-2204},{16221,-2304},{16206,-2404},{16191,-2503},{16175,-2602},{16159,-2701},{16142,-2801},{16125,-2900},{16107,-2998},{16088,-3097},{16069,-3196},{16049,-3294},{16028,-3393},{16007,-3491},{15985,-3589},{15963,-3687},{15940,-3785},{15917,-3883},{15892,-3980},{15868,-4078},{15842,-4175},{15817,-4272},{15790,-4369},{15763,-4466},{15735,-4563},{15707,-4659},{15678,-4756},{15649,-4852},{15618,-4948},{15588,-5043},{15557,-5139},{15525,-5234},{15492,-5329},{15459,-5424},{15426,-5519},{15392,-5614},{15357,-5708},{15322,-5802},{15286,-5896},{15249,-5990},{15212,-6083},{15175,-6176},{15136,-6269},{15098,-6362},{15058,-6455},{15018,-6547},{14978,-6639},{14937,-6731},{14895,-6822},{14853,-6914},{14810,-7005},{14767,-7095},{14723,-7186},{14679,-7276},{14634,-7366},{14589,-7456},{14543,-7545},{14496,-7634},{14449,-7723},{14401,-7811},{14353,-7900},{14304,-7988},{14255,-8075},{14205,-8162},{14155,-8249},{14104,-8336},{14053,-8423},{14001,-8509},{13948,-8594},{13895,-8680},{13842,-8765},{13788,-8850},{13733,-8934},{13678,-9018},{13622,-9102},{13566,-9185},{13510,-9268},{13452,-9351},{13395,-9434},{13337,-9516},{13278,-9597},{13219,-9679},{13159,-9759},{13099,-9840},{13038,-9920},{12977,-10000},{12916,-10079},{12854,-10159},{12791,-10237},{12728,-10315},{12665,-10393},{12600,-10471},{12536,-10548},{12471,-10625},{12406,-10701},{12340,-10777},{12273,-10853},{12207,-10928},{12139,-11002},{12072,-11077},{12003,-11150},{11935,-11224},{11866,-11297},{11796,-11370},{11726,-11442},{11656,-11513},{11585,-11585},{11513,-11656},{11442,-11726},{11370,-11796},{11297,-11866},{11224,-11935},{11150,-12003},{11077,-12072},{11002,-12139},{10928,-12207},{10853,-12273},{10777,-12340},{10701,-12406},{10625,-12471},{10548,-12536},{10471,-12600},{10393,-12665},{10315,-12728},{10237,-12791},{10159,-12854},{10079,-12916},{10000,-12977},{9920,-13038},{9840,-13099},{9759,-13159},{9679,-13219},{9597,-13278},{9516,-13337},{9434,-13395},{9351,-13452},{9268,-13510},{9185,-13566},{9102,-13622},{9018,-13678},{8934,-13733},{8850,-13788},{8765,-13842},{8680,-13895},{8594,-13948},{8509,-14001},{8423,-14053},{8336,-14104},{8249,-14155},{8162,-14205},{8075,-14255},{7988,-14304},{7900,-14353},{7811,-14401},{7723,-14449},{7634,-14496},{7545,-14543},{7456,-14589},{7366,-14634},{7276,-14679},{7186,-14723},{7095,-14767},{7005,-14810},{6914,-14853},{6822,-14895},{6731,-14937},{6639,-14978},{6547,-15018},{6455,-15058},{6362,-15098},{6269,-15136},{6176,-15175},{6083,-15212},{5990,-15249},{5896,-15286},{5802,-15322},{5708,-15357},{5614,-15392},{5519,-15426},{5424,-15459},{5329,-15492},{5234,-15525},{5139,-15557},{5043,-15588},{4948,-15618},{4852,-15649},{4756,-15678},{4659,-15707},{4563,-15735},{4466,-15763},{4369,-15790},{4272,-15817},{4175,-15842},{4078,-15868},{3980,-15892},{3883,-15917},{3785,-15940},{3687,-15963},{3589,-15985},{3491,-16007},{3393,-16028},{3294,-16049},{3196,-16069},{3097,-16088},{2998,-16107},{2900,-16125},{2801,-16142},{2701,-16159},{2602,-16175},{2503,-16191},{2404,-16206},{2304,-16221},{2204,-16234},{2105,-16248},{2005,-16260},{1905,-16272},{1805,-16284},{1705,-16294},{1605,-16305},{1505,-16314},{1405,-16323},{1305,-16331},{1205,-16339},{1105,-16346},{1004,-16353},{904,-16359},{803,-16364},{703,-16368},{603,-16372},{502,-16376},{402,-16379},{301,-16381},{201,-16382},{100,-16383},{0,-16384},{-100,-16383},{-201,-16382},{-301,-16381},{-402,-16379},{-502,-16376},{-603,-16372},{-703,-16368},{-803,-16364},{-904,-16359},{-1004,-16353},{-1105,-16346},{-1205,-16339},{-1305,-16331},{-1405,-16323},{-1505,-16314},{-1605,-16305},{-1705,-16294},{-1805,-16284},{-1905,-16272},{-2005,-16260},{-2105,-16248},{-2204,-16234},{-2304,-16221},{-2404,-16206},{-2503,-16191},{-2602,-16175},{-2701,-16159},{-2801,-16142},{-2900,-16125},{-2998,-16107},{-3097,-16088},{-3196,-16069},{-3294,-16049},{-3393,-16028},{-3491,-16007},{-3589,-15985},{-3687,-15963},{-3785,-15940},{-3883,-15917},{-3980,-15892},{-4078,-15868},{-4175,-15842},{-4272,-15817},{-4369,-15790},{-4466,-15763},{-4563,-15735},{-4659,-15707},{-4756,-15678},{-4852,-15649},{-4948,-15618},{-5043,-15588},{-5139,-15557},{-5234,-15525},{-5329,-15492},{-5424,-15459},{-5519,-15426},{-5614,-15392},{-5708,-15357},{-5802,-15322},{-5896,-15286},{-5990,-15249},{-6083,-15212},{-6176,-15175},{-6269,-15136},{-6362,-15098},{-6455,-15058},{-6547,-15018},{-6639,-14978},{-6731,-14937},{-6822,-14895},{-6914,-14853},{-7005,-14810},{-7095,-14767},{-7186,-14723},{-7276,-14679},{-7366,-14634},{-7456,-14589},{-7545,-14543},{-7634,-14496},{-7723,-14449},{-7811,-14401},{-7900,-14353},{-7988,-14304},{-8075,-14255},{-8162,-14205},{-8249,-14155},{-8336,-14104},{-8423,-14053},{-8509,-14001},{-8594,-13948},{-8680,-13895},{-8765,-13842},{-8850,-13788},{-8934,-13733},{-9018,-13678},{-9102,-13622},{-9185,-13566},{-9268,-13510},{-9351,-13452},{-9434,-13395},{-9516,-13337},{-9597,-13278},{-9679,-13219},{-9759,-13159},{-9840,-13099},{-9920,-13038},{-10000,-12977},{-10079,-12916},{-10159,-12854},{-10237,-12791},{-10315,-12728},{-10393,-12665},{-10471,-12600},{-10548,-12536},{-10625,-12471},{-10701,-12406},{-10777,-12340},{-10853,-12273},{-10928,-12207},{-11002,-12139},{-11077,-12072},{-11150,-12003},{-11224,-11935},{-11297,-11866},{-11370,-11796},{-11442,-11726},{-11513,-11656},{-11585,-11585},{-11656,-11513},{-11726,-11442},{-11796,-11370},{-11866,-11297},{-11935,-11224},{-12003,-11150},{-12072,-11077},{-12139,-11002},{-12207,-10928},{-12273,-10853},{-12340,-10777},{-12406,-10701},{-12471,-10625},{-12536,-10548},{-12600,-10471},{-12665,-10393},{-12728,-10315},{-12791,-10237},{-12854,-10159},{-12916,-10079},{-12977,-10000},{-13038,-9920},{-13099,-9840},{-13159,-9759},{-13219,-9679},{-13278,-9597},{-13337,-9516},{-13395,-9434},{-13452,-9351},{-13510,-9268},{-13566,-9185},{-13622,-9102},{-13678,-9018},{-13733,-8934},{-13788,-8850},{-13842,-8765},{-13895,-8680},{-13948,-8594},{-14001,-8509},{-14053,-8423},{-14104,-8336},{-14155,-8249},{-14205,-8162},{-14255,-8075},{-14304,-7988},{-14353,-7900},{-14401,-7811},{-14449,-7723},{-14496,-7634},{-14543,-7545},{-14589,-7456},{-14634,-7366},{-14679,-7276},{-14723,-7186},{-14767,-7095},{-14810,-7005},{-14853,-6914},{-14895,-6822},{-14937,-6731},{-14978,-6639},{-15018,-6547},{-15058,-6455},{-15098,-6362},{-15136,-6269},{-15175,-6176},{-15212,-6083},{-15249,-5990},{-15286,-5896},{-15322,-5802},{-15357,-5708},{-15392,-5614},{-15426,-5519},{-15459,-5424},{-15492,-5329},{-15525,-5234},{-15557,-5139},{-15588,-5043},{-15618,-4948},{-15649,-4852},{-15678,-4756},{-15707,-4659},{-15735,-4563},{-15763,-4466},{-15790,-4369},{-15817,-4272},{-15842,-4175},{-15868,-4078},{-15892,-3980},{-15917,-3883},{-15940,-3785},{-15963,-3687},{-15985,-3589},{-16007,-3491},{-16028,-3393},{-16049,-3294},{-16069,-3196},{-16088,-3097},{-16107,-2998},{-16125,-2900},{-16142,-2801},{-16159,-2701},{-16175,-2602},{-16191,-2503},{-16206,-2404},{-16221,-2304},{-16234,-2204},{-16248,-2105},{-16260,-2005},{-16272,-1905},{-16284,-1805},{-16294,-1705},{-16305,-1605},{-16314,-1505},{-16323,-1405},{-16331,-1305},{-16339,-1205},{-16346,-1105},{-16353,-1004},{-16359,-904},{-16364,-803},{-16368,-703},{-16372,-603},{-16376,-502},{-16379,-402},{-16381,-301},{-16382,-201},{-16383,-100},};

void radix2_dit(input_window<cint16> * x_in,output_window<cint16> * y_out);