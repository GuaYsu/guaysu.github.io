---
title: "从零理解大地电磁测深｜循序渐进的自学教程"
date: 2026-10-21
permalink: /posts/2026/10/mt-course/
excerpt: "从场和三角函数起步，逐步理解复数、微分方程、正演、观测、信号处理、张量、反演与地质解释。每一个新数学符号都有来历。"
tags: [大地电磁, 地球物理, 自学教程, 数学基础]
---
{% include mt-course-math.html %}

## 为什么要写这一套教程？

大地电磁测深（MT）从地表的电磁场记录推断地下导电结构，串联起电磁学、数学、信号处理、数值计算和地质学。它最难的地方之一，是许多知识的前后联系十分紧密：如果还不知道“偏导数是什么”，就要理解“对电场取旋度”，接下来的公式几乎会全部失去意义。

这套教程的目标是让**数学基础暂时薄弱的人，也能循序渐进学到完整的理论与应用**。进度不必快；有些原本放在一页里的推导，我们宁愿用整整一节解释。

每篇尽可能遵循这样的节奏：从具体的物理问题出发，明确所有新符号，再用非常简单的数字练习数学，逐步引入一般符号，写出严格公式，最后说明物理意义与适用条件。末尾的自检题用于发现自己是否真的弄懂，而非检验记忆公式的速度。

![从天然电磁场观测到地下电性结构与地质解释的流程](/images/posts/mt-course/06-data-workflow.svg)

## 第一阶段：让抽象物理量变得可见

1. **[场究竟是什么：从风场到地表电磁观测](/posts/2026/10/mt-course-01-field-wave/)**：标量、矢量、场的空间与时间依赖，$E_x$ 下标的意义，场和波的区别。
2. **[三角函数里的振幅、频率与相位](/posts/2026/10/mt-course-01b-trigonometry/)**：弧度、周期、Hz、角频率以及相位超前和滞后，带完整数字例题。
3. **[复数从哪里来：把振幅和相位装在一起](/posts/2026/10/mt-course-01c-complex/)**：$i^2=-1$、复平面、模长、辐角、欧拉公式、复振幅和复指数求导。

先把这三篇弄清楚，再进入真正的微分方程。不要跳过后面的计算题：能看懂讲解和能够自己做一次，是两件事。

## 第二阶段：研究电磁场为什么会这样变化

4. **[导数与偏导数：从坡度开始](/posts/2026/10/mt-course-02a-derivatives/)**：逐步计算变化率，理解二阶导数。
5. **[旋度与麦克斯韦方程](/posts/2026/10/mt-course-02b-curl/)**：从小回路开始理解空间旋度。
6. **[一维控制方程与指数解](/posts/2026/10/mt-course-02c-1d-ode/)**：消元、两次求导、指数函数试解。
7. **[阻抗、45° 相位与趋肤深度](/posts/2026/10/mt-course-02d-impedance/)**：从电磁场求比值。
8. **[层状地球的阻抗递推](/posts/2026/10/mt-course-02e-layers/)**：界面条件、双曲正切与自下而上递推。
9. **[有限差分与二维、三维数值正演专题](/posts/2026/10/mt-course-03-finite-difference/)**：已经拆为七节，依次讲网格、三点差分、非均匀离散、矩阵、边界、二维/三维与灵敏度。
   - [连续地下为什么要切成网格？](/posts/2026/10/mt-course-03a-grid/)
   - [三点差分与泰勒展开](/posts/2026/10/mt-course-03b-differences/)
   - [非均匀网格与电磁方程离散](/posts/2026/10/mt-course-03c-nonuniform/)
   - [矩阵与线性方程组](/posts/2026/10/mt-course-03d-matrix/)
   - [边界条件与数值误差](/posts/2026/10/mt-course-03e-boundaries/)
   - [二维、三维有限差分](/posts/2026/10/mt-course-03f-multidimensional/)
   - [薄层敏感性与时域电磁场](/posts/2026/10/mt-course-03g-sensitivity/)

## 第三阶段：从真实观测中获得可靠响应

10. **[野外观测专题](/posts/2026/10/mt-course-04-field/)**：已拆为两篇，从电压表与五通道入手，再讲测站设计和质量控制。
   - [五个通道究竟在测什么？](/posts/2026/10/mt-course-04a-measure/)
   - [从研究目标到测站布设](/posts/2026/10/mt-course-04b-station/)
11. **[数据处理专题](/posts/2026/10/mt-course-05-processing/)**：已拆为五篇，从正弦波叠加走向复杂的稳健阻抗估计。
   - [傅里叶分解与频率投影](/posts/2026/10/mt-course-05a-fourier/)
   - [采样、混叠与 FFT](/posts/2026/10/mt-course-05b-sampling/)
   - [切窗、去趋势与频谱泄漏](/posts/2026/10/mt-course-05c-windowing/)
   - [最小二乘阻抗估计](/posts/2026/10/mt-course-05d-least-squares/)
   - [远参考与稳健回归](/posts/2026/10/mt-course-05e-robust/)
12. **[阻抗张量与数据分析专题](/posts/2026/10/mt-course-06-analysis/)**：已拆成八篇，从矩阵和旋转开始，深入维度、畸变、相位张量与静态位移。
   - [阻抗张量和 tipper](/posts/2026/10/mt-course-06a-tensor/)
   - [1D、2D、3D 与 TE/TM](/posts/2026/10/mt-course-06b-dimensionality/)
   - [坐标旋转与不变量](/posts/2026/10/mt-course-06c-rotation/)
   - [WAL 不变量与 skew](/posts/2026/10/mt-course-06d-invariants/)
   - [GB 畸变分解](/posts/2026/10/mt-course-06e-distortion/)
   - [相位张量与椭圆](/posts/2026/10/mt-course-06f-phase-tensor/)
   - [走向分析方法](/posts/2026/10/mt-course-06g-strike/)
   - [静态位移与 MTpy](/posts/2026/10/mt-course-06h-static-shift/)

## 第四阶段：推断地下结构并解释地质成因

13. **[地球物理反演与 Occam 专题](/posts/2026/10/mt-course-07-inversion/)**：已整理为八篇循序渐进的讲义，从最简单的反问题与零空间开始，逐步推导 Gauss–Newton 和 Occam，最后学习反演验证。
   - [反演为什么比正演难？](/posts/2026/10/mt-course-07a-inverse-problem/)
   - [零空间和非唯一性](/posts/2026/10/mt-course-07b-nullspace/)
   - [误差、残差与归一化 RMS](/posts/2026/10/mt-course-07c-misfit/)
   - [雅可比矩阵与灵敏度](/posts/2026/10/mt-course-07d-jacobian/)
   - [梯度、Hessian 与 Newton](/posts/2026/10/mt-course-07e-newton/)
   - [Gauss–Newton 逐步推导](/posts/2026/10/mt-course-07f-gauss-newton/)
   - [Occam、模型粗糙度与 L 曲线](/posts/2026/10/mt-course-07g-occam/)
   - [分辨率、棋盘测试与模型可靠性](/posts/2026/10/mt-course-07h-validation/)
14. **[电性结构与地质解释专题](/posts/2026/10/mt-course-08-interpretation/)**：已拆成九篇，从岩石物性出发，逐步分析大陆、俯冲带、古裂谷、SinoProbe、西准噶尔和海洋 MT 案例，最后建立完整证据链。
   - [电阻率、电子与离子传导](/posts/2026/10/mt-course-08a-conduction/)
   - [Archie 定律与流体连通性](/posts/2026/10/mt-course-08b-archie/)
   - [大陆阵列与全球地幔电性](/posts/2026/10/mt-course-08c-global-arrays/)
   - [Cascadia 俯冲带导电异常](/posts/2026/10/mt-course-08d-cascadia/)
   - [北美裂谷、缝合带及大陆深部成像](/posts/2026/10/mt-course-08e-northamerica/)
   - [SinoProbe：华北和青藏高原](/posts/2026/10/mt-course-08f-sinoprobe/)
   - [西准噶尔 MT 与高密度 AMT](/posts/2026/10/mt-course-08g-junggar-amt/)
   - [海洋 MT：海沟、洋中脊与熔体](/posts/2026/10/mt-course-08h-marine/)
   - [如何写有证据的地质解释](/posts/2026/10/mt-course-08i-evidence/)

## 统一的坐标与时间约定

为了避免不同章节中出现“同一符号代表不同方向”的混乱，本系列采用如下**基准坐标约定**：

- $x$ 轴正向指**北**，$y$ 轴正向指**东**，$z$ 轴正向指**下**（深度增大方向）；这三个轴组成右手坐标系（NED）。
- $E_x$ 是电场的北向分量，$E_y$ 是东向分量。电极实际铺设方向可以有所偏离，但必须记录方位并按约定旋转。
- 一维正演里的深度 $z\ge 0$ 向地下增加；传播的衰减解写为 $e^{-kz}$，选择 $\operatorname{Re}k>0$。
- 单频相量使用 $e^{+i\omega t}$ 时间因子，因此 $\partial/\partial t\mapsto+i\omega$；法拉第定律写作 $\nabla\times\mathbf E=-i\omega\mu_0\mathbf H$。
- 讨论二维构造时，常为了方便把局部的 $x'$ 轴旋转到地质走向，此时带撇号的方向不必仍指地理正北。旋转采用本系列约定的 $\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T$。

以上只是明确的计算约定，不是声称全世界的资料都这样设置。有些处理软件采用东—北—上（ENU）或反号时间因子，使用外部数据或公式前应核查坐标、单位与傅里叶约定。

## 统一的数学记号

- $E,H$ 为电场和磁场强度，$E_x,H_y$ 表示在相应坐标方向的分量。
- $t$ 为时间，$f$ 为频率，$T$ 为周期，$\omega=2\pi f$ 为角频率。
- $i^2=-1$，$\widetilde E$ 表示复振幅，$\operatorname{Re}$ 表示取实部。
- $\rho$ 为电阻率，$\sigma=1/\rho$ 为电导率；$\mu_0$ 为真空磁导率。
- 当研究单个频率的场时，默认使用 $e^{i\omega t}$ 时间约定；变化时会明确标注。

不同章节出现的新符号仍会在首次使用时重新解释，避免要求读者把这里当作一张必须背下的符号表。

## 延伸专题

[从解析解到数值解：为什么“算出来了”还不够](/posts/2026/09/analytical-vs-numerical-solution/) 详细讨论一维正演中的数值误差、网格加密、边界和收敛检验。它更适合在学完有限差分之后阅读。

> **修订说明**：本系列已按照零基础、循序渐进的路线完成所有主要主题的细化，从场与数学基础直到电性结构及地质解释。各单元保留原有专题链接，便于站内引用与后续修订。