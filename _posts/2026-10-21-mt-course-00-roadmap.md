---
title: "大地电磁测深系统讲义｜导读与学习路线"
date: 2026-10-21
permalink: /posts/2026/10/mt-course/
excerpt: "从电磁场的直觉开始，循序学习数学、MT 正演、野外观测、信号处理、张量分析、反演和地质解释；附八份课件索引。"
tags: [大地电磁, 地球物理, 学习讲义, MT]
---
{% include mt-course-math.html %}

> **怎样使用这套讲义？** 从第 1 讲往下读；遇到公式先看每个字母表示什么，再看假设、推导与物理意义，最后动手算例。内容依据 Bo YANG 老师 2026 年《Fundamentals of Geophysics》电磁学部分八份课件重新组织，增加数学解释和教学图示。课件与本系列之间不是逐页对应；每节末尾标注重要页码。

![课程全流程](/images/posts/mt-course/06-data-workflow.svg)

## 目录

1. [第 1 讲：场、波、复数与趋肤深度](/posts/2026/10/mt-course-01-field-wave/)
2. [第 2 讲：麦克斯韦方程到一维解析正演](/posts/2026/10/mt-course-02-maxwell-forward/)
3. [第 3 讲：有限差分与二维、三维正演](/posts/2026/10/mt-course-03-finite-difference/)
4. [第 4 讲：野外仪器、观测设计与数据质量](/posts/2026/10/mt-course-04-field/)
5. [第 5 讲：时间序列、FFT、阻抗估计与稳健处理](/posts/2026/10/mt-course-05-processing/)
6. [第 6 讲：维度分析、张量旋转、畸变和静态位移](/posts/2026/10/mt-course-06-analysis/)
7. [第 7 讲：反演、最小二乘、雅可比与 Occam](/posts/2026/10/mt-course-07-inversion/)
8. [第 8 讲：从电阻率结构走向地质解释](/posts/2026/10/mt-course-08-interpretation/)

## 核心问题与八份课件的对应

| 课件 | 本系列覆盖重点 |
|---|---|
| Ch.1 Introduction（50页） | MT 概念、天然场源、方法对比、工作流程 |
| Ch.2 Theory Part I（35页） | 麦克斯韦方程、1D 控制方程、阻抗、层状介质 |
| Ch.2 Theory Part II（32页） | FDM、1D/2D/3D 正演、敏感性、时域波动 |
| Ch.3 MT in the Wild（33页） | 仪器、选点、布设、噪声与野外安全 |
| Ch.4 Data Processing（55页） | 切窗、FFT、采样、TF、远参考与稳健回归 |
| Ch.5 Data Analysis（68页） | 维度、旋转、WAL、GB、相位张量、静态位移、MTpy |
| Ch.6 Data Inversion（69页） | 反演非唯一性、最小二乘、GN、Occam、软件与诊断 |
| Ch.7 Interpretation（110页） | 导电机制、全球尺度、大陆 MT 阵列、案例、海洋 MT |

每章采用“直觉 → 数学工具 → 课件公式 → 物理解释 → 实际应用”的次序。**一维、二维、三维的‘维数’描述的是地下电性结构随几个空间方向改变，不是电磁场有几个分量。**这是贯穿全课的重要概念。

## 学习前约定

- 复数电磁场采用课件中的时间因子 $e^{+i\omega t}$，所以 $\partial/\partial t\to i\omega$。不同书籍采用相反的时间约定，公式中虚部符号可能相反。
- 地下深度 $z$ 向下为正。$E_x$ 与 $H_y$ 表示电、磁场在各坐标轴上的分量，而非“电场与磁场只有两种方向”。
- 电导率 $\sigma$ 单位 S/m，电阻率 $\rho=1/\sigma$ 单位 $\Omega\cdot$m。角频率 $\omega=2\pi f$。
- 先理解理想化模型，再说明何时不再适用；真实地质解释必须结合地质、地震、岩石物性与分辨率信息。

## 进一步阅读

站内已有专题：[从解析解到数值解：为什么“算出来了”还不够](/posts/2026/09/analytical-vs-numerical-solution/)。该文深入讨论网格、边界、地表导数与数值误差；本讲义先补齐理解它的数学与物理基础。