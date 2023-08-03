# CCFSys定制计算挑战赛2023：决赛

> 参赛学校：华东师范大学
> 
> 指导教师：徐飞老师
> 
> 参赛队员：白卓岩、堵仪萱、徐珑珊

## 工作概述

本分支是大赛**决赛设计**的提交。在初赛作品的基础上，决赛阶段优化了1K、4K、8K点数的AIE设计以及8K点数PL和AIE的连接。

本项目实现的是FFT赛题，从初赛到决赛的整体完成情况如下。

- 完成了1024点`cint16`、`cfloat`类型的小规模FFT算法的AIE设计、优化和仿真；
- 探索了AIE kernel之间的连接，扩展到更大点数（4k、8k）的`cint16`类型的FFT；
- 实现了PL和AIE的数据连接，可以从host端调用，能够在VCK5000硬件上运行。

## 相关文档
- 答辩PPT：[答辩PPT](https://github.com/abuqiqi/CCC2023/blob/main/%E7%AD%94%E8%BE%A9PPT.pptx)。
- 设计报告：[设计报告](https://github.com/abuqiqi/CCC2023/blob/main/%E8%AE%BE%E8%AE%A1%E6%8A%A5%E5%91%8A.pdf)。
- 初赛提交：见本仓库的[preliminary分支](https://github.com/abuqiqi/CCC2023/tree/preliminary)。
- 工作日志：见本仓库的[backup分支](https://github.com/abuqiqi/CCC2023/tree/backup)。

## 代码运行
决赛阶段的代码是本仓库`sources`文件夹下提交的`fft_8k`工程，可以执行AIE仿真和硬件运行。

1. AIE仿真

执行以下代码可以单独运行AIE的仿真。
> 注：执行`make`之前需要确认`Makefile`文件中`XPFM`路径。

```shell
# 编译并运行AIE仿真
cd ./CCC2023/sources/fft_8k/aie
make
make aieemu
```

2. 硬件运行

在`sources/fft_8k/execution`文件夹下存放了通过主机调用PL和AIE必要的`fft.xclbin`文件、`host.exe`文件和输入文件`DataInFFTO.txt`，以及运行完毕所产生的输出文件`DataOutFFT0.txt`。如需在VCK5000上运行，可执行以下代码。

```shell
# 克隆hacc_demo仓库
git clone https://github.com/Xtra-Computing/hacc_demo.git

# 获取VCK5000计算节点（根据hacc_demo存放路径修改指令）
./hacc_demo/env/vck5000_alloc 3
source ./hacc_demo/env/vck5000_env

# 在本项目的execution文件夹下运行可执行文件
./CCC2023/sources/fft_4k/execution/host.exe

# 退出节点
./hacc_demo/env/vck5000_exit
```

执行完毕后，可使用`sources/fft_8k/notebook`文件夹下的`.ipynb`文件可视化输出结果并进行验证。

## 目录说明
决赛提交的主要目录结构如下。
```
CCC2023
├── sources
│   ├── fft_8k          8K-point FFT完整代码
│   │   ├── aie         AIE设计
│   │   ├── execution   可执行文件
│   │   ├── host        host端代码
│   │   ├── hw_link     硬件连接
│   │   ├── notebook    数据生成、可视化
│   │   ├── pl          PL设计
│   │   └── Makefile
│   │
│   ├── fft_1k          1K-point FFT AIE代码
│   └── fft_4k          4K-point FFT AIE代码
│
├── README.md
├── 答辩PPT.pptx
└── 设计报告.pdf
```