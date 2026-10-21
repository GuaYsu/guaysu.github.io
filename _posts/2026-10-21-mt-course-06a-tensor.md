---
title: "数据分析 01｜阻抗张量为什么是四个复数：从矩阵乘法开始"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-06a-tensor/
excerpt: "从一个方向的比例关系到两个磁场输入、两个电场输出，逐项理解阻抗矩阵和垂直磁场传递函数。"
tags: [大地电磁, 阻抗张量, 数学基础]
---
{% include mt-course-math.html %}

前面我们已经用时间序列、傅里叶变换和最小二乘得到了每个频率的电场与磁场复数系数。但有个问题：为什么一个理想一维地球只需 $Z_{xy}=E_x/H_y$，一到真实三维地下，就出现四个 $Z$？

## 一、先从一个输入和一个输出开始

设一个频率上的复数磁场是 $H_y$，电场是 $E_x$，理想模型满足

$$E_x=Z_{xy}H_y.$$

这只是“输出=比例系数×输入”。例如 $H_y=2$、$Z_{xy}=3+4i$，则 $E_x=6+8i$。复数系数同时包含振幅比例和相位关系。

现在地下可能沿水平方向变化。$E_x$ 不一定只由 $H_y$ 控制，也可能与 $H_x$ 有关。因此用两个系数：

$$E_x=Z_{xx}H_x+Z_{xy}H_y.$$

同理，另一个输出 $E_y$ 也需要两个系数：

$$E_y=Z_{yx}H_x+Z_{yy}H_y.$$

两行叠在一起，就成了阻抗张量。

![从两路磁场到两路电场](/images/posts/mt-course/21-impedance-matrix.svg)

## 二、矩阵只是两行方程的简写

$$
\boxed{
\begin{pmatrix}E_x\\E_y\end{pmatrix}=
\begin{pmatrix}Z_{xx}&Z_{xy}\\Z_{yx}&Z_{yy}\end{pmatrix}
\begin{pmatrix}H_x\\H_y\end{pmatrix}.
}
$$

矩阵乘法按“每一行乘以输入列，然后相加”计算。第一行结果就是 $Z_{xx}H_x+Z_{xy}H_y$。下标第一位标记**输出电场方向**，第二位标记**输入磁场方向**。这些下标从来都不是测站位置。

## 三、给四个系数随便设数，算一遍

先拿实数练习：

$$
\mathbf Z=\begin{pmatrix}1&2\\3&4\end{pmatrix},
\quad
\mathbf H=\begin{pmatrix}5\\6\end{pmatrix}.
$$

第一行 $E_x=1\times5+2\times6=17$；第二行 $E_y=3\times5+4\times6=39$，所以

$$
\mathbf E=\begin{pmatrix}17\\39\end{pmatrix}.
$$

真实阻抗的四项一般是复数，但计算规则完全一样。**一个 $2\times2$ 复矩阵含有四个复数，也就是八个实数**。之后的维度分析，就是研究八个数之间是否存在特殊关系。

## 四、磁场用 $H$ 还是 $B$？

如果数据文件里记录的是磁感应强度 $\mathbf B$，且采用 $\mathbf B=\mu_0\mathbf H$，那么将 $\mathbf E=\mathbf Z_H\mathbf H$ 改写为

$$
\mathbf E=\underbrace{\frac{\mathbf Z_H}{\mu_0}}_{\mathbf M}\mathbf B.
$$

因此用 $E/H$ 和 $E/B$ 定义的传递系数数值及单位不同。处理软件和教材若采用不同定义，必须先换算后比较。

## 五、第五个通道 $H_z$ 如何使用？

除了水平的四个阻抗元素，还可以估计**垂直磁场传递函数**：

$$
H_z=T_xH_x+T_yH_y.
$$

$T_x,T_y$ 组成 tipper（倾子、垂直传递函数）向量。它反映电磁感应响应的方向性和某些横向异常；它不是把地下电阻率“直接量出来”的工具。绘制 induction arrows（感应箭头）时，还要确认所用箭头约定，有些约定方向相反。

## 六、先认识特殊的一维模型

对水平层状、各向同性且场源可以用平面波近似的模型：

$$
\mathbf Z=
\begin{pmatrix}0&Z\\-Z&0\end{pmatrix}.
$$

这意味着 $Z_{xx}=Z_{yy}=0$，且 $Z_{xy}=-Z_{yx}$。两个对角元消失不是由于“测错了”，而是介质水平方向的对称性与电磁场的极化结构决定的。

**自检：** 如果 $\mathbf Z=\begin{pmatrix}1&2\\3&4\end{pmatrix}$、$\mathbf H=(0,5)^T$，输出是多少？答：$(10,20)^T$。记住只有第二列对结果起作用。

[下一篇：一维、二维、三维与 TE/TM](/posts/2026/10/mt-course-06b-dimensionality/) · [专题目录](/posts/2026/10/mt-course-06-analysis/)