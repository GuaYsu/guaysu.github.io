---
title: "数据分析 03｜为什么旋转阻抗张量，要左右各乘一个矩阵？"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-06c-rotation/
excerpt: "从二维箭头在不同坐标系的投影讲起，推导旋转矩阵、转置和 Z'=RZR^T，附具体数字算例。"
tags: [大地电磁, 线性代数, 坐标旋转, 阻抗张量]
---
{% include mt-course-math.html %}

终于轮到一个看起来很吓人的公式：

$$
\boxed{\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T}.
$$

我们要做的不是背它，而是回答：**既然空间中的电、磁场根本没变，只是换了一套坐标，为什么左右会各出现一个旋转矩阵？**

## 一、向量本身不变，分量可以改变

![同一个电场在不同坐标系中的分量](/images/posts/mt-course/22-coordinate-rotation.svg)

设一根箭头指向东北方，长度为 $\sqrt2$。若 $x$ 指北、$y$ 指东，它有分量

$$\mathbf v=\begin{pmatrix}1\\1\end{pmatrix}.$$

现在把新 $x'$ 轴转到箭头方向，新的 $y'$ 轴与其垂直；箭头当然还是原来的箭头，但它在新坐标里的分量会成为 $(\sqrt2,0)$。

所以**坐标变换改变的是分量数字，不是物理场本身**。

## 二、旋转矩阵从两个投影得到

图中的横纵轴用于说明抽象的二维分量投影，不对应纸面上的地理东、西方向。全系列仍采用 $x$ 北、$y$ 东、$z$ 下的右手坐标系。设新坐标轴相对旧坐标旋转角度 $\theta$，选定如下被动坐标变换约定：

$$
\mathbf R(\theta)=
\begin{pmatrix}
\cos\theta&\sin\theta\\
-\sin\theta&\cos\theta
\end{pmatrix}.
$$

那么

$$
\begin{pmatrix}v_x'\\v_y'\end{pmatrix}
=
\mathbf R
\begin{pmatrix}v_x\\v_y\end{pmatrix},
$$

展开得

$$
v_x'=v_x\cos\theta+v_y\sin\theta,\qquad
v_y'=-v_x\sin\theta+v_y\cos\theta.
$$

$45^\circ$ 时 $\cos\theta=\sin\theta=1/\sqrt2$，把 $(1,1)$ 代入，就得到 $v_x'=\sqrt2,v_y'=0$。我们用具体数字验证了旋转公式。

## 三、什么叫转置？为什么旋转的逆矩阵特别简单？

矩阵转置写作上标 $T$：把行和列互换。例如

$$
\begin{pmatrix}a&b\\c&d\end{pmatrix}^T
=\begin{pmatrix}a&c\\b&d\end{pmatrix}.
$$

对旋转矩阵，

$$
\mathbf R^T=
\begin{pmatrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{pmatrix}.
$$

一个有趣的性质是 $\mathbf R^T\mathbf R=\mathbf I$，这里 $\mathbf I=\begin{pmatrix}1&0\\0&1\end{pmatrix}$ 叫单位矩阵，乘上它不会改变向量。

为什么？第一行乘第一列得到 $\cos^2\theta+\sin^2\theta=1$；第一行乘第二列得到 $-\cos\theta\sin\theta+\sin\theta\cos\theta=0$，其余两项同理。所以

$$\mathbf R^{-1}=\mathbf R^T.$$

这句话的直觉：坐标转过去以后，再转回来，结果应与原来完全相同。

## 四、真正推导阻抗张量旋转公式

旧坐标中的物理关系是

$$\mathbf E=\mathbf Z\mathbf H.$$

现在同一个电场和磁场换成新坐标：

$$
\mathbf E'=\mathbf R\mathbf E,\qquad
\mathbf H'=\mathbf R\mathbf H.
$$

由第二式得到 $\mathbf H=\mathbf R^{-1}\mathbf H'=\mathbf R^T\mathbf H'$。把它代回电场关系：

$$
\begin{aligned}
\mathbf E'
&=\mathbf R\mathbf E\\
&=\mathbf R\mathbf Z\mathbf H\\
&=\mathbf R\mathbf Z\mathbf R^T\mathbf H'.
\end{aligned}
$$

新坐标中的阻抗定义为 $\mathbf E'=\mathbf Z'\mathbf H'$。比较两边，得到

$$
\boxed{\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T}.
$$

**左边的 $\mathbf R$ 把输出电场转到新坐标；右边的 $\mathbf R^T$ 把新磁场输入还原到旧坐标，以便使用旧坐标的阻抗关系。** 不是随便在两侧各加一个矩阵。

## 五、算一个 90° 的具体例子

设

$$
\mathbf Z=\begin{pmatrix}a&b\\c&d\end{pmatrix},\qquad
\mathbf R_{90}=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
$$

先计算

$$
\mathbf R_{90}\mathbf Z
=\begin{pmatrix}c&d\\-a&-b\end{pmatrix}.
$$

再右乘

$$
\mathbf R_{90}^T=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
$$

得到

$$
\boxed{\mathbf Z'=\begin{pmatrix}d&-c\\-b&a\end{pmatrix}}.
$$

如果原来只有非对角元，新坐标中仍只有非对角元，但它们的位置和符号会变化。以后看到文献里 $Z_{xy}$ 与 $Z_{yx}$ 符号不同，一定要先确认坐标与时间约定。

## 六、最后认识三个旋转不变量

有些数字虽然矩阵分量变了，却始终保持不变。对于水平正交旋转，最常用的有：

$$
\operatorname{tr}\mathbf Z=Z_{xx}+Z_{yy}\quad\text{（迹）},
$$

$$
\operatorname{atr}\mathbf Z=Z_{xy}-Z_{yx}\quad\text{（反迹）},
$$

$$
\det\mathbf Z=Z_{xx}Z_{yy}-Z_{xy}Z_{yx}\quad\text{（行列式）}.
$$

为什么迹不变？因为 $\operatorname{tr}(\mathbf R\mathbf Z\mathbf R^T)=\operatorname{tr}(\mathbf Z\mathbf R^T\mathbf R)=\operatorname{tr}\mathbf Z$，用到了矩阵迹的循环置换性质。行列式不变则因为 $\det(\mathbf R)=1$，所以 $\det\mathbf Z'=\det\mathbf R\det\mathbf Z\det\mathbf R^T=\det\mathbf Z$。

反迹的旋转不变性可以利用二维反对称矩阵检验；它对上述保持方向的旋转不变，但若坐标变换包含镜像反射，要特别处理方向符号。

## 七、自检

1. $\mathbf v=(1,0)^T$ 转动坐标 $90^\circ$ 后是什么？答：$(0,-1)^T$，在本节约定下。
2. 为什么 $\mathbf Z'$ 不是简单的 $\mathbf R\mathbf Z$？因为磁场输入也要变换坐标。
3. 迹与反迹相同是否足以证明两个阻抗张量物理相同？不够，它们只包含部分信息。

[下一节：WAL 不变量与维度判断](/posts/2026/10/mt-course-06d-invariants/) · [专题目录](/posts/2026/10/mt-course-06-analysis/)