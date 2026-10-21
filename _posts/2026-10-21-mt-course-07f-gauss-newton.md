---
title: "地球物理反演 06｜Gauss–Newton 从哪里来：每个矩阵都推给你看"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-07f-gauss-newton/
excerpt: "从非线性最小二乘和一阶泰勒展开推导 Gauss–Newton 正规方程，手算 F(m)=m² 的两次迭代，再加入模型正则化。"
tags: [大地电磁, Gauss Newton, 反演算法, 数学推导]
---
{% include mt-course-math.html %}

到这里，我们已经学会了数据残差、雅可比矩阵、梯度和 Hessian。终于可以解答一个核心问题：

**为什么 Gauss–Newton（高斯–牛顿）方法不需要把复杂的正演函数完整求二阶导数，却仍能进行 Newton 型更新？**

我们从最简单的最小二乘目标开始，每一步都写出来。

## 一、先把数据误差写成一个向量

令观测数据为 $\mathbf d$，预测数据为 $F(\mathbf m)$。定义

$$
\mathbf r(\mathbf m)=F(\mathbf m)-\mathbf d.
$$

这是残差向量：预测减观测。为简化推导，先不加数据权重，目标函数写作

$$
\boxed{\Phi_d(\mathbf m)=\frac12\mathbf r^T\mathbf r.}
$$

为什么前面多了 $1/2$？因为求导时平方项会产生系数 2，乘上 $1/2$ 后恰好消掉。它只改变目标函数的整体数值尺度，不改变最小值位置。

## 二、用一阶泰勒展开预测下一步的残差

当前模型为 $\mathbf m_k$，想改变 $\Delta\mathbf m$。一阶展开：

$$
F(\mathbf m_k+\Delta\mathbf m)
\approx F(\mathbf m_k)+\mathbf J_k\Delta\mathbf m.
$$

其中 $\mathbf J_k$ 是当前模型的雅可比矩阵。

减去观测：

$$
\mathbf r_{\rm new}
\approx\mathbf r_k+\mathbf J_k\Delta\mathbf m.
$$

下一步的近似目标函数因此变为

$$
\Phi_{\rm approx}(\Delta\mathbf m)
=\frac12\|\mathbf r_k+\mathbf J_k\Delta\mathbf m\|^2.
$$

这一步是整个 Gauss–Newton 的关键：我们没有对原本很复杂的非线性正演反复求二阶导数，而是把它**在当前模型附近暂时线性化**。

## 三、先把向量平方展开

向量平方范数表示向量各元素平方的总和。利用矩阵乘法：

$$
\begin{aligned}
2\Phi_{\rm approx}
&=(\mathbf r+\mathbf J\Delta\mathbf m)^T
(\mathbf r+\mathbf J\Delta\mathbf m)\\
&=\mathbf r^T\mathbf r
+2\Delta\mathbf m^T\mathbf J^T\mathbf r
+\Delta\mathbf m^T\mathbf J^T\mathbf J\Delta\mathbf m.
\end{aligned}
$$

第二行的三个部分分别是：原来的残差平方、随模型改动出现的一次项、模型改动产生的二次项。就像普通代数中的

$$(a+b)^2=a^2+2ab+b^2.$$

区别只是这里的 $a,b$ 是向量，乘法要保留矩阵的顺序。

## 四、对模型增量求导

对上面的二次函数对 $\Delta\mathbf m$ 求导，得到

$$
\frac{\partial\Phi_{\rm approx}}{\partial(\Delta\mathbf m)}
=\mathbf J^T\mathbf r+\mathbf J^T\mathbf J\Delta\mathbf m.
$$

令其为零：

$$
\mathbf J^T\mathbf J\Delta\mathbf m
=-\mathbf J^T\mathbf r.
$$

因此

$$
\boxed{
\mathbf J^T\mathbf J\Delta\mathbf m=-\mathbf J^T\mathbf r.
}
$$

这就是无权重、无正则化情况下的 Gauss–Newton 正规方程。左边的 $\mathbf J^T\mathbf J$ 是一个**近似 Hessian**。实际求解时不一定计算它的逆矩阵，可以用 QR、SVD 或适当的迭代线性求解器。

## 五、亲手做两次 Gauss–Newton 更新

设只有一个参数，正演函数是

$$
F(m)=m^2,
$$

实际观测为 $d=9$。我们想寻找令 $m^2\approx9$ 的模型。先从 $m_0=2$ 出发。

此时预测 $F(2)=4$，残差

$$r_0=4-9=-5.$$

雅可比在一元情况下就是导数：

$$J_0=\frac{dF}{dm}\bigg|_{m=2}=2m=4.$$

正规方程变成

$$J_0^2\Delta m=-J_0r_0.$$

代入：

$$16\Delta m=-4(-5)=20,$$

所以

$$\Delta m=1.25,\qquad m_1=2+1.25=3.25.$$

新的预测为

$$F(3.25)=10.5625.$$

它已经比旧预测的 4 更接近观测 9。现在重新计算雅可比：

$$J_1=2(3.25)=6.5.$$

更新量：

$$
\Delta m=-\frac{10.5625-9}{6.5}\approx-0.24038.
$$

所以

$$
m_2=3.25-0.24038\approx3.00962.
$$

我们正一步步接近正解 $m=3$。同时这个例子也提醒我们：若从 $m_0=0$ 出发，$J_0=0$，最小二乘正规方程便失去有效的更新方向。初始模型会影响求解。

## 六、加入数据权重后怎样变化？

实际数据误差不等，因此令 $\mathbf W_d$ 是数据白化或加权矩阵，使

$$
\Phi_d=\frac12\|\mathbf W_d\mathbf r\|^2.
$$

按照同样的展开步骤，可以得到

$$
\boxed{
\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf J\,\Delta\mathbf m
=-\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf r.
}
$$

这里 $\mathbf W_d^T\mathbf W_d$ 对应数据协方差逆矩阵（在适当正定条件下）。它的作用是让误差较小的数据更有效地约束模型，同时正确处理相关误差。

## 七、再加入模型正则化

为防止参数沿弱约束方向剧烈震荡，可以使用

$$
\Phi=
\frac12\|\mathbf W_d(F(\mathbf m)-\mathbf d)\|^2
+\frac{\lambda}{2}\|\mathbf L(\mathbf m-\mathbf m_{\rm ref})\|^2.
$$

$\lambda$ 为正则化权重，$\mathbf L$ 是模型约束算子，$\mathbf m_{\rm ref}$ 是参考模型。

在当前模型附近线性化、求导、令零，得到

$$
\boxed{
\begin{aligned}
(&\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf J
+\lambda\mathbf L^T\mathbf L)\Delta\mathbf m
={}&-\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf r\\
&-\lambda\mathbf L^T\mathbf L(\mathbf m_k-\mathbf m_{\rm ref}).
\end{aligned}
}
$$

这个式子很长，但你已经知道它每部分来自哪里：

- $\mathbf J$ 来自正演函数的一阶导；
- $\mathbf r$ 来自预测与观测的差；
- $\mathbf W_d$ 来自观测误差；
- $\mathbf L$ 来自对模型形态的要求；
- $\lambda$ 决定我们对拟合与模型约束的权衡。

解出增量以后，再使用步长 $\gamma$ 更新：

$$
\boxed{\mathbf m_{k+1}=\mathbf m_k+\gamma\Delta\mathbf m.}
$$

这就是地球物理反演常见的“正演 → 误差 → 灵敏度 → 更新 → 新正演”迭代框架。

## 八、它为什么不直接使用完整 Hessian？

对无权重的最小二乘目标，精确 Hessian 含有

$$
\mathbf H_{\rm exact}
=\mathbf J^T\mathbf J
+\sum_j r_j\,\nabla^2F_j.
$$

第二项涉及预测函数的**二阶导数**。Gauss–Newton 近似忽略这项，使用 $\mathbf J^T\mathbf J$。当残差较小等条件满足时，这往往是合理而计算上更可行的做法；但在强非线性、大残差问题里也可能遇到困难。

下一节再看一个更关键的问题：既然正则化已经出现，Occam 反演为何会特别强调“在满足数据拟合的模型中尽量平滑”？

[继续：Occam 与正则化的完整来历](/posts/2026/10/mt-course-07g-occam/) · [反演专题目录](/posts/2026/10/mt-course-07-inversion/)