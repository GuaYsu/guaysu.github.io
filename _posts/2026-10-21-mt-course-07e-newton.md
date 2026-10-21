---
title: "地球物理反演 05｜梯度和 Hessian 在干什么：从一元二次函数学 Newton 方法"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-07e-newton/
excerpt: "用一个抛物线手算导数与二阶导数，再推导 Newton 更新，理解梯度、Hessian 和步长搜索。"
tags: [大地电磁, 梯度, Hessian, Newton方法]
---
{% include mt-course-math.html %}

反演通常需要让某个误差指标尽量小。我们已经会计算失配，也懂得“改变地下参数会怎样改变预测数据”。现在只差最后一步：

**知道某个模型的误差比较大，计算机如何决定下一个模型往哪边改？**

我们从最简单的抛物线开始。

## 一、假设目标函数是一条抛物线

设

$$
\Phi(m)=(m-3)^2.
$$

$m$ 是我们要调整的唯一参数。展开为

$$
\Phi(m)=m^2-6m+9.
$$

这是开口向上的抛物线。它在 $m=3$ 时恰好等于零，其余位置均为正，因此最小值就在 $m=3$。

但假设我们不知道这个结论，只能利用导数寻找最低点，应该怎么做？

## 二、一阶导数告诉我们往哪边下降

求导：

$$
\frac{d\Phi}{dm}=2(m-3).
$$

如果当前 $m=0$，导数为 $-6$，意味着沿 $m$ 增大的方向，$\Phi$ 正在降低。反过来如果 $m=5$，导数为 $+4$，则向左走可以降低目标函数。

当

$$
\frac{d\Phi}{dm}=0
$$

时，$m=3$。这时目标函数附近的斜率为零，它可能是局部最低点、最高点或其他驻点，还需要进一步判断。

## 三、二阶导数告诉我们曲线弯曲得怎样

再对一阶导数求导：

$$
\frac{d^2\Phi}{dm^2}=2.
$$

二阶导数大于 0，说明曲线局部朝上弯；在一元情形下，如果一阶导数为零、二阶导数大于零，则它是局部极小点。

对多个模型参数，

$$
\mathbf m=\begin{pmatrix}m_1\\m_2\\\vdots\\m_M\end{pmatrix},
$$

所有一阶偏导组成**梯度向量**：

$$
\mathbf g=\nabla\Phi=
\begin{pmatrix}
\partial\Phi/\partial m_1\\
\partial\Phi/\partial m_2\\
\vdots\\
\partial\Phi/\partial m_M
\end{pmatrix}.
$$

二阶偏导组成**Hessian 矩阵**：

$$
\boxed{H_{jk}=\frac{\partial^2\Phi}{\partial m_j\partial m_k}.}
$$

不要把这里的粗体 Hessian $\mathbf H$ 与之前的磁场强度 $\mathbf H$ 混淆。它们恰好用了相同字母，具体含义要看上下文。

## 四、Newton 方法为什么需要二阶信息？

在当前模型 $\mathbf m$ 附近，把目标函数展开到二阶：

$$
\Phi(\mathbf m+\Delta\mathbf m)
\approx
\Phi(\mathbf m)
+\mathbf g^T\Delta\mathbf m
+\frac12\Delta\mathbf m^T\mathbf H\Delta\mathbf m.
$$

这行公式其实只是在把“当前位置的函数值、斜率和弯曲程度”加在一起。对一维情况，它就是熟悉的

$$
\Phi(m+\Delta m)\approx
\Phi(m)+\Phi'(m)\Delta m+\frac12\Phi''(m)(\Delta m)^2.
$$

把上面的二阶近似对 $\Delta\mathbf m$ 求导，并令结果为零：

$$
\mathbf g+\mathbf H\Delta\mathbf m=\mathbf0.
$$

移项：

$$
\boxed{\mathbf H\Delta\mathbf m=-\mathbf g.}
$$

这就是 Newton 更新方程。它的意思是：**用当前位置的斜率和局部弯曲程度，预测怎样修改参数才能降低目标函数。**

## 五、先手算一个 Newton 更新

仍用

$$
\Phi(m)=(m-3)^2.
$$

取初始 $m_0=0$。此时

$$
g=-6,\qquad H=2.
$$

Newton 方程：

$$
2\Delta m=-(-6)=6,
$$

所以

$$
\Delta m=3,\qquad m_1=m_0+\Delta m=3.
$$

一次就到了最低点！这是因为我们的目标函数恰好是二次函数，二阶展开没有遗漏更高阶项。

实际 MT 反演的正演函数更复杂，局部二次近似只在有限范围内有效，因此通常不会一步完成。

## 六、为什么还要有步长搜索？

![沿着目标函数进行逐次更新的示意](/images/posts/mt-course/31-newton-path.svg)

即使算得一个候选增量 $\Delta\mathbf m$，也不一定立刻迈完整一步。实际模型可以按

$$
\boxed{\mathbf m_{k+1}=\mathbf m_k+\gamma\Delta\mathbf m}
$$

更新，其中 $\gamma$ 是**步长**。可以尝试 $\gamma=1,1/2,1/4,\ldots$，比较实际目标函数是否下降。

这叫**线搜索（line search）**。因为真实目标函数可能有多个局部极小值，局部方向也可能失准，步长控制可以减少因一步太大而失败的情况。

Newton 型方法也不天然保证收敛到全局最佳解。初始模型、参数尺度、正则化、目标函数形状都会影响迭代。

## 七、两道练习

**①** $\Phi(m)=(m-5)^2$，当前 $m=1$。梯度为 $-8$，二阶导数为 2，Newton 增量为 $\Delta m=4$，新模型是 5。

**②** 如果 Hessian 很接近奇异矩阵，直接求逆是否可靠？往往不可靠。它说明某些模型方向的曲率很小，可能与零空间或病态问题有关，需要正则化与稳定求解策略。

下一节我们把本节的 Newton 更新与前面的雅可比连接起来，得到地球物理反演里更常用的 Gauss–Newton 方法。

[继续：Gauss–Newton 为什么不必完整计算二阶正演导数？](/posts/2026/10/mt-course-07f-gauss-newton/) · [反演专题目录](/posts/2026/10/mt-course-07-inversion/)