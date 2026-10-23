---
title: "数据分析 05｜近地表导体怎样扭曲电场：Groom–Bailey 畸变"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-06e-distortion/
excerpt: "从导电小体旁的电流和界面电荷开始，理解电流性畸变矩阵、twist、shear、anisotropy 与 site gain。"
tags: [大地电磁, Groom Bailey, 电流性畸变, 数据分析]
---
{% include mt-course-math.html %}

想象我们真正关心的是十几千米深处的大型导电结构。可测站正好落在一小块浅部异常导体旁边，这个小体的尺度远小于区域构造，却会改变局部电流的方向和强弱。

于是仪器测到的电场，可能是**区域电场经过局部结构重新分配以后形成的结果**。

![浅部导体与局部电流畸变](/images/posts/mt-course/24-galvanic-distortion.svg)

## 一、为什么“近地表异常”可能产生电场畸变？

电导率不同的两个材料相接时，电流连续性和电磁边界条件可能导致界面上积累电荷，从而产生额外的局部电场。若异常体尺度远小于研究频带的特征电磁变化尺度，局部响应中的一部分可以近似看成频率无关的**电流性畸变**（galvanic distortion）。

在这种理想近似里：

$$
\boxed{\mathbf E_{\rm obs}=\mathbf C\mathbf E_{\rm reg}}.
$$

$\mathbf E_{\rm reg}$ 是我们希望知道的区域电场；$\mathbf E_{\rm obs}$ 是测到的电场；$\mathbf C$ 是 $2\times2$ **实数矩阵**。

举个简单例子：

$$
\mathbf C=\begin{pmatrix}2&0\\0&1\end{pmatrix}.
$$

若区域电场是 $(3,4)^T$，那么观测电场就成为 $(6,4)^T$。你看，原本相同的物理区域场，在测点上被不同方向的比例放大了。

## 二、为什么阻抗也会被畸变？

如果暂时假设磁场不发生明显畸变，且区域电场满足

$$
\mathbf E_{\rm reg}=\mathbf Z_{\rm reg}\mathbf H,
$$

代入上一式：

$$
\mathbf E_{\rm obs}
=\mathbf C\mathbf Z_{\rm reg}\mathbf H.
$$

因此测得的阻抗是

$$
\boxed{\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg}}.
$$

这就解释了一个令人困惑的现象：区域地下本来近似二维，但测到的阻抗矩阵可能有明显的对角项。原因未必是深部构造突然变成三维，也可能是浅部电场发生了方向混合。

## 三、GB 分解试图怎样认识矩阵 C？

Groom–Bailey（GB）分解把局部畸变视为若干具有较直观作用的运算组合。常见符号写作

$$
\boxed{\mathbf C=g\,\mathbf T\,\mathbf S\,\mathbf A},
$$

其中 $g$ 为 site gain（测站增益），$\mathbf T$ 为 twist（扭转），$\mathbf S$ 为 shear（剪切），$\mathbf A$ 为 anisotropy（各向异性增益）。不同文献使用的规范化和矩阵排列约定需要核对。

为了先理解数学，取一种常用的示意参数化。扭转可写作一个旋转类矩阵：

$$
\mathbf T(t)=\frac1{\sqrt{1+t^2}}
\begin{pmatrix}1&-t\\t&1\end{pmatrix}.
$$

它表示两个水平电场分量因局部效应而发生类似旋转的混合。参数 $t$ 控制这种“扭”的程度；其对应角度可由反正切表示。

剪切可以示意为

$$
\mathbf S(s)=\frac1{\sqrt{1+s^2}}
\begin{pmatrix}1&s\\s&1\end{pmatrix}.
$$

于是如果原区域电场为 $(E_x,E_y)^T$，剪切后的两个分量为

$$
E_x'=\frac{E_x+sE_y}{\sqrt{1+s^2}},\qquad
E_y'=\frac{sE_x+E_y}{\sqrt{1+s^2}}.
$$

两条式子告诉我们：**测到的每个方向，都混入了原本另一个方向的场**。

各向异性增益可以用类似

$$
\mathbf A(a)=\frac1{\sqrt{1-a^2}}
\begin{pmatrix}1+a&0\\0&1-a\end{pmatrix},\qquad |a|<1
$$

的形式表示，它使两个方向的幅值按不同的比例变化。标量 $g$ 则使整个电场统一放大或缩小。上述因子可以采用不同的归一化形式，部分比例会被吸收入 $g$。

## 四、为什么把四个效应分开仍很难唯一求出？

假设只在一个频率、一个测站测到了四个复阻抗元素，共八个实数。我们不仅想恢复区域二维阻抗的两个复数分量，还要估计区域走向、局部扭转、剪切以及增益等参数。

未知量过多，而且某些因子可互相补偿。例如整体放大 $g$ 与区域阻抗整体模长的改变，在单站数据里就可能无法区分。这意味着**测到一个完整阻抗张量，并不保证我们能唯一恢复全部畸变参数**。

常见策略是跨多个频率与多个测点联合拟合，并对局部畸变的频率无关性、区域二维走向的一致性等提出约束。即便拟合效果好，绝对 site gain 通常仍可能存在不可辨识问题。

## 五、畸变不是一种单一现象

电流性（galvanic）畸变侧重界面电荷及局部电流重新分配；电感性（inductive）畸变与频率相关的感应电流、次生磁场有关。

同一个浅部导体对不同频率和 TE/TM 模式可能产生不同影响。将所有异常响应简单叫作“静态位移”，会混淆幅值变化、方向混合与频率相关的效应。

**自检：** 若 $\mathbf C=\begin{pmatrix}1&0.5\\0&1\end{pmatrix}$，区域场为 $(2,4)^T$，观测场是多少？答：$(4,4)^T$。因为 $E_x$ 混进了 $0.5E_y$。

[下一节：相位张量为何能抵消实数电流性畸变？](/posts/2026/10/mt-course-06f-phase-tensor/) · [专题目录](/posts/2026/10/mt-course-06-analysis/)