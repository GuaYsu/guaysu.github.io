---
title: "反演入门 05｜高斯–牛顿法为什么要解一个矩阵方程？"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-07e-gauss-newton/
excerpt: "从抛物线最低点、梯度和曲率学起，一步一步推导普通牛顿法、Gauss–Newton 与正则化模型更新。"
tags: [地球物理反演, Gauss Newton, Hessian, 线搜索]
---
{% include mt-course-math.html %}

上一节，我们认识了雅可比矩阵：它告诉我们改变地下模型参数，预测数据会朝哪个方向变化。可是反演的目标是让数据失配减小，怎样利用这些变化率找到更好的模型？

答案之一是**高斯–牛顿（Gauss–Newton）方法**。名字听着很吓人，我们仍然先从一条最简单的抛物线开始。

## 一、一个数的最小值怎样找？

假设

$$
\Phi(m)=(m-3)^2.
$$

这是开口向上的抛物线，在 $m=3$ 处最小，最小值为零。我们可以用导数找最低点：

$$
\frac{d\Phi}{dm}=2(m-3).
$$

令导数等于零：

$$
2(m-3)=0\quad\Rightarrow\quad m=3.
$$

为什么要让导数为零？因为在平滑函数的内部局部极小值处，往左和往右的一阶变化倾向会平衡，切线斜率为零。

## 二、梯度是什么？

当参数从一个 $m$ 变成很多个 $m_1,m_2,\dots,m_K$ 时，我们无法只写一个普通导数，而要对每个参数分别求偏导：

$$
\mathbf g=\nabla\Phi=
\begin{pmatrix}
\partial\Phi/\partial m_1\\
\partial\Phi/\partial m_2\\
\vdots\\
\partial\Phi/\partial m_K
\end{pmatrix}.
$$

$\mathbf g$ 叫**梯度向量**，它记录每个参数稍微增加时，目标函数最初会怎样变化。它指向局部上升最快的方向；负梯度则可作为下降方向。

## 三、Hessian 又是什么？

一阶导数描述坡度，二阶导数描述坡度变化的快慢。在多参数情况下，把所有二阶偏导数排成一个矩阵：

$$
\mathbf H_{jk}=
\frac{\partial^2\Phi}{\partial m_j\partial m_k}.
$$

这个矩阵叫 Hessian，描述目标函数在当前模型附近的**曲率**。比如一个方向非常陡，目标函数对它敏感；另一个方向非常平，改变模型可能对目标函数影响不大。

利用二阶泰勒展开，目标函数在当前模型附近可近似为

$$
\Phi(\mathbf m+\Delta\mathbf m)
\approx \Phi(\mathbf m)
+\mathbf g^T\Delta\mathbf m
+\frac12\Delta\mathbf m^T\mathbf H\Delta\mathbf m.
$$

右边依次是：原目标值、一阶线性变化、二阶曲率修正。对模型增量 $\Delta\mathbf m$ 求导并令其为零，得到

$$
\boxed{\mathbf H\Delta\mathbf m=-\mathbf g}.
$$

这便是多参数牛顿更新方程的来历。它仍然只是“找近似抛物面的最低点”，只是现在变成了矩阵形式。

## 四、高斯–牛顿怎样用雅可比避免复杂的二阶正演导数？

先考虑数据失配

$$
\Phi_d=\frac12\|\mathbf W_d(F(\mathbf m)-\mathbf d)\|^2.
$$

$\mathbf d$ 为观测，$F(\mathbf m)$ 为预测，$\mathbf W_d$ 是数据误差权重矩阵。

定义残差向量 $\mathbf r=F(\mathbf m)-\mathbf d$。小幅修改模型后，借助雅可比线性化：

$$
F(\mathbf m+\Delta\mathbf m)-\mathbf d
\approx \mathbf r+\mathbf J\Delta\mathbf m.
$$

因此局部近似目标函数成为

$$
\Phi_d(\mathbf m+\Delta\mathbf m)
\approx\frac12\|\mathbf W_d(\mathbf r+\mathbf J\Delta\mathbf m)\|^2.
$$

现在问题变成了关于 $\Delta\mathbf m$ 的**二次函数最小化**。求导并令其为零，得到普通高斯–牛顿的近似法方程：

$$
\boxed{
\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf J\,\Delta\mathbf m
=-\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf r.
}
$$

与真正的 Hessian 相比，$\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf J$ 忽略了一部分与正演函数二阶导数有关的项，使计算更容易。这就是“Gauss–Newton”与完整 Newton 法的主要区别之一。

这里用的是实数数据向量形式；复数阻抗可拆为实部和虚部后使用实数形式，也可采用相应复共轭处理。

## 五、加入正则化以后，为什么左右各多出一项？

上一节已经学过模型约束：

$$
\Phi_m=\frac12\|\mathbf L(\mathbf m-\mathbf m_{\rm ref})\|^2.
$$

加上权重 $\lambda$，总目标为

$$
\Phi=\Phi_d+\lambda\Phi_m.
$$

令 $\mathbf q=\mathbf m-\mathbf m_{\rm ref}$。模型更新后，参考模型差异变为 $\mathbf q+\Delta\mathbf m$。因此局部优化问题为

$$
\frac12\|\mathbf W_d(\mathbf r+\mathbf J\Delta\mathbf m)\|^2+
\frac\lambda2\|\mathbf L(\mathbf q+\Delta\mathbf m)\|^2.
$$

分别求导，将含 $\Delta\mathbf m$ 的项移至左侧，得到

$$
\boxed{
(\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf J+
\lambda\mathbf L^T\mathbf L)\Delta\mathbf m
=-\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf r
-\lambda\mathbf L^T\mathbf L\mathbf q.
}
$$

这就是常见的**正则化高斯–牛顿更新方程**。左侧第一项与数据对模型的敏感性有关，第二项与模型平滑约束有关；右边告诉我们当前预测误差与偏离先验的情况。

不用一遍背下它。你应能够逐项说出：哪个来自数据、哪个来自模型约束、哪个是待求的模型变化。

## 六、算出增量后为什么不能直接全部加上去？

线性化是在当前模型附近近似成立的。如果更新太大，新的模型可能远离原来的局部近似区域，导致目标函数反而增大。

因此可以引入步长 $\gamma$：

$$
\boxed{\mathbf m_{\rm new}=\mathbf m+\gamma\Delta\mathbf m}.
$$

常见做法是在多个候选 $\gamma$ 上计算目标函数，选择能让拟合和约束合理改善的步长。这称为**线搜索**。

完整的一次迭代就是：当前模型正演、计算残差与雅可比、解更新方程、选步长、更新模型，再从新的模型重新开始。大型反演常需更先进的线性求解与优化策略，但基本逻辑相通。

## 七、自检

1. $\Phi(m)=(m-5)^2$ 的导数是什么？答：$2(m-5)$。
2. 雅可比与 Hessian 的差别是什么？雅可比描述预测数据对模型参数的一阶变化率；Hessian 描述目标函数对模型参数的二阶曲率。
3. 为什么需要线搜索？因为局部线性近似不能保证任意大步更新都有利。

[下一节：L 曲线、恢复试验和怎样判断反演结果可信](/posts/2026/10/mt-course-07f-diagnostics/) · [专题目录](/posts/2026/10/mt-course-07-inversion/)