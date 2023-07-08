# CCFSys定制计算挑战赛2023

> 参赛学校：华东师范大学
> 
> 指导教师：徐飞老师
> 
> 参赛队员：白卓岩、堵仪萱、徐珑珊

## 工作概述

本项目实现的是FFT赛题，完成情况如下。

- 完成了1024点`cint16`、`cfloat`类型的小规模FFT算法的AIE设计、优化和仿真；
- 探索了AIE kernel之间的连接，扩展到更大点数（4k、8k）的`cint16`类型的FFT；
- 实现了PL和AIE的数据连接，可以从host端调用，完成了系统级的仿真并能够在VCK5000硬件上运行。

## 相关文档
- 设计报告：[设计报告](https://github.com/abuqiqi/CCC2023/blob/main/%E8%AE%BE%E8%AE%A1%E6%8A%A5%E5%91%8A.pdf)。
- 工作日志：见本仓库的[backup分支](https://github.com/abuqiqi/CCC2023/tree/backup)。

## 代码运行
1. AIE仿真

执行以下代码可以单独运行AIE的仿真。

```shell
# 编译并运行AIE仿真
cd ./CCC2023/sources/fft_4k/aie
make
make aieemu
```

2. 硬件运行

在`sources/fft_4k/execution`文件夹下存放了通过主机调用PL和AIE必要的`fft.xclbin`文件、`host.exe`文件和输入文件`DataInFFTO.txt`，以及运行完毕所产生的输出文件`DataOutFFT0.txt`。如需在VCK5000上运行，可执行以下代码。

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

执行完毕后，可使用`sources/fft_4k/notebook`文件夹下的`fft.ipynb`验证输出结果。
