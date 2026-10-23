---
title: "13｜雅可比、Newton 与 Gauss–Newton 全推导｜合并完整版"
date: 2026-10-21
modified: 2026-10-23
permalink: /posts/2026/10/mt-course-07b-optimization/
excerpt: "相关章节完整合并，保留全部原文、数学过程、实例和自检问题。"
tags: [大地电磁, 数学基础, 自学教程]
---
{% include mt-course-math.html %}

## 本篇目录
- [反演入门 04｜雅可比矩阵是什么：改动一层地下，数据会改变多少？](#sec-07d-jacobian)
- [地球物理反演 05｜梯度和 Hessian 在干什么：从一元二次函数学 Newton 方法](#sec-07e-newton)
- [地球物理反演 06｜Gauss–Newton 从哪里来：每个矩阵都推给你看](#sec-07f-gauss-newton)


---

<a id="sec-07d-jacobian"></a>
## 第1单元｜反演入门 04｜雅可比矩阵是什么：改动一层地下，数据会改变多少？

前面我们已经有了目标：寻找能够解释观测、同时不过于复杂的地下电阻率模型。那么在计算机里，模型应该怎么修改？如果随便把某一层电阻率变大、另一层变小，很容易改错方向。

有一个特别重要的工具，可以告诉我们：**修改模型中的某个参数，会影响哪些数据、影响多大。** 它叫雅可比矩阵（Jacobian matrix），也叫灵敏度矩阵。

### 一、从一个简单偏导数开始

先不看电磁学，考虑函数

$$
d_1=2m_1+3m_2.
$$

这里 $m_1,m_2$ 是两个可调参数；$d_1$ 是预测数据。若只增加 $m_1$，保持 $m_2$ 不变，$d_1$ 每增加一个单位就增加 2，所以

$$
\frac{\partial d_1}{\partial m_1}=2.
$$

同理

$$
\frac{\partial d_1}{\partial m_2}=3.
$$

回忆偏导数的意思：**改动一个参数时，暂时固定其他参数，看看结果变化有多快。**

现在再增加两条预测：

$$
d_2=m_1-m_2,\qquad d_3=4m_2.
$$

求偏导数后，可以排成表：

| 数据 | 对 $m_1$ 的偏导 | 对 $m_2$ 的偏导 |
|---|---:|---:|
| $d_1$ | 2 | 3 |
| $d_2$ | 1 | −1 |
| $d_3$ | 0 | 4 |

这个表就是雅可比矩阵：

$$
\boxed{\mathbf J=
\begin{pmatrix}
2&3\\
1&-1\\
0&4
\end{pmatrix}}.
$$

![雅可比矩阵的行和列](/images/posts/mt-course/29-jacobian.svg)

注意它的行列意义：**一行对应一条数据，一列对应一个模型参数**。这里有三条数据、两个模型参数，所以 $\mathbf J$ 有三行、两列。

### 二、改动参数以后，怎样预估数据改变？

设模型参数的改变量为

$$
\Delta\mathbf m=
\begin{pmatrix}\Delta m_1\\\Delta m_2\end{pmatrix}.
$$

如果 $\Delta m_1=0.1$、$\Delta m_2=-0.2$，就把它们放进刚才的表：

$$
\Delta\mathbf d=\mathbf J\Delta\mathbf m.
$$

逐行做乘法：

$$
\Delta d_1=2(0.1)+3(-0.2)=-0.4,
$$

$$
\Delta d_2=1(0.1)-1(-0.2)=0.3,
$$

$$
\Delta d_3=0(0.1)+4(-0.2)=-0.8.
$$

因此

$$
\Delta\mathbf d=
\begin{pmatrix}-0.4\\0.3\\-0.8\end{pmatrix}.
$$

我们把两个模型参数的改动，转换成了三条预测数据的改动。这就是雅可比矩阵最直观的工作。

### 三、真实 MT 正演不是简单加法，怎么办？

对于真实地层，电阻率与阻抗关系通常是非线性的。预测数据写作

$$
\mathbf d_{\rm pred}=F(\mathbf m).
$$

这里 $F$ 是完整的正演程序。即使无法写出简单的显式函数，我们仍然可以在某个模型 $\mathbf m$ 附近，研究小扰动怎样改变预测。

对第 $j$ 条数据和第 $k$ 个参数定义

$$
\boxed{J_{jk}=\frac{\partial F_j(\mathbf m)}{\partial m_k}}.
$$

它的物理含义与刚才完全相同：固定其余参数，只把第 $k$ 个参数微调一下，观察第 $j$ 条数据对它有多敏感。

### 四、从高中学过的切线过渡到线性化

对于一个变量，我们学过

$$
f(x+\Delta x)\approx f(x)+f'(x)\Delta x.
$$

这个式子是在说：**在当前位置附近，用切线预测小幅移动后函数的值。**

多个参数、多个数据时，就把一个导数换成整个雅可比矩阵：

$$
\boxed{
F(\mathbf m+\Delta\mathbf m)
\approx F(\mathbf m)+\mathbf J\Delta\mathbf m.
}
$$

右边第一项是当前模型的预测；第二项表示模型改动预计引起的响应变化。这叫**一阶泰勒线性化**。

非常重要的限制：当 $\Delta\mathbf m$ 很大、非线性很强时，这个近似可能明显失效。因此反演程序往往要多次更新模型，每次重新正演并重新计算或更新灵敏度。

### 五、如何让计算机近似求偏导数？

最朴素的方法就是**有限差分扰动**。先运行一次正演，得到 $F(\mathbf m)$。然后只把第 $k$ 个参数增加很小的 $\epsilon$，再次正演：

$$
J_{jk}\approx
\frac{
F_j(\mathbf m+\epsilon\mathbf e_k)-F_j(\mathbf m)
}{\epsilon}.
$$

$\mathbf e_k$ 是第 $k$ 个位置为 1、其余位置为 0 的向量，所以 $\epsilon\mathbf e_k$ 只改变第 $k$ 个参数。

计算办法与两点差值近似导数相同，只是现在改变的变量换成了地下模型参数。

这种方法容易理解、便于检验，但对每个参数都重新正演可能非常昂贵。大型二维、三维反演经常使用伴随法等更高效的灵敏度计算技术。

### 六、雅可比大就一定能识别参数吗？

不能。若某列数值很小，说明当前数据对相应参数变化不敏感，通常很难约束它；但若两列几乎成比例，即使各自不小，仍难区分两个参数的独立作用。

例如两列都是 $(1,2,3)^T$ 与 $(2,4,6)^T$，改变第一个参数和第二个参数可能产生同方向的数据变化。观测只知道二者的组合，仍然存在非唯一性。

较高灵敏度意味着参数变化能够引起可观的数据响应。参数能否被单独分辨，还取决于其他参数是否产生相近响应，以及噪声和正则化的影响。

### 七、自检

1. $d=5m_1+2m_2$，对 $m_1$ 的偏导是多少？答：5。
2. 雅可比有 100 行、20 列，表示多少条数据和多少个参数？答：100 条数据、20 个参数。
3. 为什么要反复更新雅可比？因为真实正演通常非线性，模型改变以后局部导数也可能改变。

[下一节：从梯度走到高斯–牛顿更新](/posts/2026/10/mt-course-07b-optimization/#sec-07e-newton) · [返回专题目录](/posts/2026/10/mt-course-07a-inverse-theory/#sec-07-inversion)


---

<a id="sec-07e-newton"></a>
## 第2单元｜地球物理反演 05｜梯度和 Hessian 在干什么：从一元二次函数学 Newton 方法

反演通常需要让某个误差指标尽量小。我们已经会计算失配，也懂得“改变地下参数会怎样改变预测数据”。现在只差最后一步：

**知道某个模型的误差比较大，计算机如何决定下一个模型往哪边改？**

我们从最简单的抛物线开始。

### 一、假设目标函数是一条抛物线

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

### 二、一阶导数告诉我们往哪边下降

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

### 三、二阶导数告诉我们曲线弯曲得怎样

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

这里粗体 $\mathbf H$ 表示 Hessian 矩阵，前面电磁学的 $\mathbf H$ 则是磁场强度。字母 $H$ 同时兼两份差事，读公式时要看上下文。

### 四、Newton 方法为什么需要二阶信息？

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

### 五、先手算一个 Newton 更新

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

### 六、为什么还要有步长搜索？

![沿着目标函数进行逐次更新的示意](/images/posts/mt-course/31-newton-path.svg)

即使算得一个候选增量 $\Delta\mathbf m$，也不一定立刻迈完整一步。实际模型可以按

$$
\boxed{\mathbf m_{k+1}=\mathbf m_k+\gamma\Delta\mathbf m}
$$

更新，其中 $\gamma$ 是**步长**。可以尝试 $\gamma=1,1/2,1/4,\ldots$，比较实际目标函数是否下降。

这叫**线搜索（line search）**。因为真实目标函数可能有多个局部极小值，局部方向也可能失准，步长控制可以减少因一步太大而失败的情况。

Newton 型方法也不天然保证收敛到全局最佳解。初始模型、参数尺度、正则化、目标函数形状都会影响迭代。

### 七、两道练习

**①** $\Phi(m)=(m-5)^2$，当前 $m=1$。梯度为 $-8$，二阶导数为 2，Newton 增量为 $\Delta m=4$，新模型是 5。

**②** 如果 Hessian 很接近奇异矩阵，直接求逆是否可靠？往往不可靠。它说明某些模型方向的曲率很小，可能与零空间或病态问题有关，需要正则化与稳定求解策略。

下一节我们把本节的 Newton 更新与前面的雅可比连接起来，得到地球物理反演里更常用的 Gauss–Newton 方法。

[继续：Gauss–Newton 为什么不必完整计算二阶正演导数？](/posts/2026/10/mt-course-07b-optimization/#sec-07f-gauss-newton) · [反演专题目录](/posts/2026/10/mt-course-07a-inverse-theory/#sec-07-inversion)


---

<a id="sec-07f-gauss-newton"></a>
## 第3单元｜地球物理反演 06｜Gauss–Newton 从哪里来：每个矩阵都推给你看

到这里，我们已经学会了数据残差、雅可比矩阵、梯度和 Hessian。终于可以解答一个核心问题：

**为什么 Gauss–Newton（高斯–牛顿）方法不需要把复杂的正演函数完整求二阶导数，却仍能进行 Newton 型更新？**

我们从最简单的最小二乘目标开始，每一步都写出来。

### 一、先把数据误差写成一个向量

令观测数据为 $\mathbf d$，预测数据为 $F(\mathbf m)$。定义

$$
\mathbf r(\mathbf m)=F(\mathbf m)-\mathbf d.
$$

这是残差向量：预测减观测。为简化推导，先不加数据权重，目标函数写作

$$
\boxed{\Phi_d(\mathbf m)=\frac12\mathbf r^T\mathbf r.}
$$

为什么前面多了 $1/2$？因为求导时平方项会产生系数 2，乘上 $1/2$ 后恰好消掉。它只改变目标函数的整体数值尺度，不改变最小值位置。

### 二、用一阶泰勒展开预测下一步的残差

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

### 三、先把向量平方展开

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

### 四、对模型增量求导

为了看清求导过程，先观察第 $k$ 个参数增量。平方残差的导数可以写成 $\sum_j J_{jk}(r_j+\sum_\ell J_{j\ell}\Delta m_\ell)$。这正是矩阵表达式第 $k$ 个位置的值。把所有参数的导数排成一列，得到

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

### 五、亲手做两次 Gauss–Newton 更新

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

### 六、加入数据权重后怎样变化？

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

### 七、再加入模型正则化

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

### 八、它为什么不直接使用完整 Hessian？

对无权重的最小二乘目标，精确 Hessian 含有

$$
\mathbf H_{\rm exact}
=\mathbf J^T\mathbf J
+\sum_j r_j\,\nabla^2F_j.
$$

第二项涉及预测函数的**二阶导数**。Gauss–Newton 近似忽略这项，使用 $\mathbf J^T\mathbf J$。当残差较小等条件满足时，这往往是合理而计算上更可行的做法；但在强非线性、大残差问题里也可能遇到困难。

下一节再看一个更关键的问题：既然正则化已经出现，Occam 反演为何会特别强调“在满足数据拟合的模型中尽量平滑”？

[继续：Occam 与正则化的完整来历](/posts/2026/10/mt-course-07c-occam-resolution/#sec-07g-occam) · [反演专题目录](/posts/2026/10/mt-course-07a-inverse-theory/#sec-07-inversion)
