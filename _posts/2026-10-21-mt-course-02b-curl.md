---
title: "大地电磁入门 05｜旋度并不可怕：从一个小回路走到麦克斯韦方程"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-02b-curl/
excerpt: "先从四条边的环流和二维小例子理解旋度，再辨认梯度、散度和麦克斯韦方程；逐一解释每一个符号。"
tags: [大地电磁, 数学基础, 旋度, 麦克斯韦方程]
---
{% include mt-course-math.html %}

在上一节，我们学会了 $\partial E_x/\partial z$：**电场的 $x$ 方向分量，随深度 $z$ 改变得有多快。**

现在电磁学要问一个略难的问题：电场在空间各个方向都有分量，这些分量又可能随着位置变化；有没有一个数学工具可以把这些信息组合起来，描述场在小回路附近的环流性质？

答案是**旋度**。先不要被名字吓到。我们从一条小正方形路径开始。

## 一、沿着一个小方框走一圈

想象空气中某处有一个速度场，你拿着小风速计，沿正方形四条边走一圈。

我们关注沿每条边方向的风速分量。比如某条边长 1 米，风沿这条边方向吹的速度是 3 个单位，那么这条边对“沿路径的流动”贡献可以记作 3；如果方向相反，贡献就记为负数。

沿整个闭合路径，把四条边的这种贡献加起来，便得到一种叫作**环流**的量。对于连续的矢量场，更严格的写法是

$$
\oint_C\mathbf E\cdot d\boldsymbol\ell.
$$

这行符号逐个看：

- $\oint_C$：沿闭合曲线 $C$ 积分，也就是绕一圈累加；
- $\mathbf E$：当前位置的电场矢量；
- $d\boldsymbol\ell$：沿路径的一小段位移矢量；
- 中间的点 $\cdot$：点积，提取电场在路径方向上的分量。

整圈的环流还不是旋度。**旋度对应把这个小圈缩得无限小之后，单位面积的有向环流密度**，并且它是一个矢量。

## 二、试试看：只有 $E_x$，为什么仍可能有旋度？

设电场只有 $x$ 方向分量，且它的大小会随深度 $z$ 变化：

$$
\mathbf E=(E_x(z),0,0).
$$

比如深度为 1 米的位置 $E_x=1$，深度为 2 米处 $E_x=2$，方向都沿正 $x$。

![电场沿 x 方向，但强度随深度改变，导致小回路的环流](/images/posts/mt-course/08-curl-loop.svg)

在 $x$–$z$ 平面里画一个小方框。如果沿正 $x$ 走上边，再沿正 $z$ 走右边，接着沿负 $x$ 走下边，最后沿负 $z$ 返回：

- 上边贡献约为 $E_x(z_{\rm top})\Delta x$；
- 两条竖边因为电场没有 $z$ 分量，贡献为 0；
- 下边反向走，其贡献约为 $-E_x(z_{\rm bottom})\Delta x$。

所以总环流近似为

$$
\big[E_x(z_{\rm top})-E_x(z_{\rm bottom})\big]\Delta x.
$$

而小方框的面积是 $\Delta x\,\Delta z$。将环流除以面积：

$$
\frac{E_x(z_{\rm top})-E_x(z_{\rm bottom})}{\Delta z}.
$$

当 $\Delta z$ 越来越小时，这个表达式趋向

$$
-\frac{\partial E_x}{\partial z}.
$$

这正是旋度在 $y$ 方向的分量之一。这里的正负号依赖于小回路的绕行方向；我们选的路线从正 $y$ 方向观察对应固定的正向约定。

所以**电场箭头本身不必真的绕圈排列，旋度仍然可能不为零**。只要相邻位置的场在小回路上无法完全相互抵消，就可能出现非零的环流。

## 三、完整旋度公式长什么样？

对三维矢量场 $\mathbf E=(E_x,E_y,E_z)$，旋度记作

$$
\nabla\times\mathbf E.
$$

其中 $\nabla$（nabla）是一个用来组织空间微分的符号，$\times$ 这里表示旋度运算，不是两个普通数字相乘。

完整展开为

$$
\nabla\times\mathbf E=
\begin{pmatrix}
\dfrac{\partial E_z}{\partial y}-\dfrac{\partial E_y}{\partial z}\\[6pt]
\dfrac{\partial E_x}{\partial z}-\dfrac{\partial E_z}{\partial x}\\[6pt]
\dfrac{\partial E_y}{\partial x}-\dfrac{\partial E_x}{\partial y}
\end{pmatrix}.
$$

对我们刚才的特例，$E_y=E_z=0$，而 $E_x$ 只依赖 $z$，所以

$$
\nabla\times\mathbf E=
\begin{pmatrix}
0\\
\dfrac{\partial E_x}{\partial z}\\
0
\end{pmatrix}.
$$

请留意：上面的完整公式规定了标准右手坐标系与正向环路方向。之前从特定绕行方向算出的有向环流可能是反号；旋度的 $y$ 分量按照右手定则定义为 $+\partial E_x/\partial z$。**不能把任意绕行顺序的环流符号，直接当作旋度分量的符号。**

还有两种常见操作，先认识即可：梯度 $\nabla f$ 描述标量场增长最快的方向；散度 $\nabla\cdot\mathbf E$ 描述局部源汇性质。与旋度不同，散度的结果是标量。在直角坐标系中

$$
\nabla\cdot\mathbf E=
\frac{\partial E_x}{\partial x}
+\frac{\partial E_y}{\partial y}
+\frac{\partial E_z}{\partial z}.
$$

它把各方向分量沿**自身方向**的变化率相加，而旋度组合的是交叉的空间变化率。

## 四、终于可以读懂第一条麦克斯韦方程

最常见的一条写作

$$
\boxed{\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t}}.
$$

我们现在能读出两边了。

左边是电场的旋度，描述局部环流性质；右边是磁感应强度 $\mathbf B$ 随时间变化的负值。

这条方程可以用一句话记住：**磁场的时间变化，会与电场的空间环流联系起来**。它是法拉第电磁感应定律的微分表达形式。

这里 $\mathbf B$ 与上一节提到的 $\mathbf H$ 不同，但在简单的线性、各向同性材料中常满足

$$
\mathbf B=\mu\mathbf H,
$$

$\mu$ 是材料的磁导率。对通常近似为非磁性的岩石，我们常取 $\mu\approx\mu_0$。不能看到 $B$ 和 $H$ 就认为是两种完全无关的东西，也不能把它们的数值单位混为一谈。

还有一条对导电地下尤其重要：

$$
\boxed{\nabla\times\mathbf H=\mathbf J+\frac{\partial\mathbf D}{\partial t}}.
$$

这里 $\mathbf J$ 是电流密度，$\mathbf D$ 是电位移矢量。$\mathbf J$ 表示导电电荷实际运动造成的电流；$\partial\mathbf D/\partial t$ 是**位移电流项**。

在某些低频、导电介质中，位移电流相对于传导电流很小，可以作近似

$$
\nabla\times\mathbf H\approx\mathbf J=\sigma\mathbf E.
$$

$\sigma$ 读作“西格玛”，是电导率，表示介质导电的容易程度。电阻率是 $\rho=1/\sigma$，它们互为倒数。这个简化将在后续推导中使用，但不能在任何频率、任何介质里无条件成立。

## 五、先用一维几何把三维式子缩小

假设地层水平方向均匀，电场选择 $x$ 方向分量，磁场选择 $y$ 方向分量，而且它们都只随深度 $z$ 改变：

$$
\mathbf E=(E_x(z),0,0),\qquad
\mathbf H=(0,H_y(z),0).
$$

由刚才的旋度公式，

$$
\nabla\times\mathbf E=
\begin{pmatrix}0\\dE_x/dz\\0\end{pmatrix},
$$

而

$$
\nabla\times\mathbf H=
\begin{pmatrix}-dH_y/dz\\0\\0\end{pmatrix}.
$$

把它们与麦克斯韦方程比较，再采用 $e^{i\omega t}$ 的简谐时间表示：

$$
\boxed{\frac{dE_x}{dz}=-i\omega\mu_0H_y},
$$

$$
\boxed{-\frac{dH_y}{dz}=\sigma E_x}.
$$

现在这两条重要的一维方程，不会再像突然冒出来的天书了。它们是三维旋度方程，在非常明确的模型条件下，抽取对应分量之后的结果。

注意这里出现的 $i$ 还是复数单位；因为采用 $e^{i\omega t}$，时间偏导数在相量表示下变成了乘以 $i\omega$。有关原因我们已经在前一篇文章里推导过。

## 六、做两道最简单的理解题

**问题一：** $E_x(z)=3z$，$E_y=E_z=0$。求旋度在 $y$ 方向的分量。

根据 $(\nabla\times\mathbf E)_y=\partial E_x/\partial z-\partial E_z/\partial x$，第一项为 3，第二项为 0，结果为 3。

**问题二：** 为什么从 $\nabla\times\mathbf H=\sigma\mathbf E$ 的 $x$ 分量中会得到一个负号？

把旋度的 $x$ 分量展开：

$$
(\nabla\times\mathbf H)_x=
\frac{\partial H_z}{\partial y}
-\frac{\partial H_y}{\partial z}.
$$

在一维假设中 $H_z=0$，所以剩下 $-dH_y/dz$，因此对应 $\sigma E_x$。

如果你已经能够回答这两题，我们就具备下一步的条件：把两条一阶方程结合起来，得到电场随深度的二阶方程。

[继续：一维控制方程为什么会出现指数解](/posts/2026/10/mt-course-02c-1d-ode/) · [课程目录](/posts/2026/10/mt-course/)