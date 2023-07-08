# CCFSys定制计算挑战赛

> 参赛学校：华东师范大学
> 
> 指导教师：徐飞老师
> 
> 参赛队员：白卓岩、堵仪萱、徐珑珊

## 工作日志

- AIE仿真
    - fft_v0_dsplib：调用dsplib库函数实现1024-cint16-FFT。
    - fft_v1_cfloat_aie：自己设计算法，实现1024-cfloat-FFT。
    - fft_v2_cfloat_aie：优化算法，实现1024-cfloat-FFT。
    - fft_v2_cint16_aie：将数据精度降为cint16，使用移位处理omega，实现1024-cint16-FFT。

- AIE+PL硬件运行
    - fft_1k_pl：在fft_v2_cint16_aie的基础上，增加pl、host，可以在vck5000上运行。
    - fft_4k_v1_pl：在fft_1k_pl的基础上，修改AIE、PL，将点数扩展到4096。
    - fft_4k_v2_pl：在fft_4k_v1_pl的基础上，优化了AIE设计。