---
title: "MT 讲义 07｜从反演非唯一性到高斯–牛顿与 Occam"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-07-inversion/
excerpt: "理解正演与反演的差别，建立最小二乘目标函数，解释雅可比、梯度、Hessian、正则化、Occam 与数据拟合权衡。"
tags: [大地电磁, 地球物理反演, 最小二乘, Occam]
---
{% include mt-course-math.html %}

## 循序渐进的六篇反演数学教程

1. [01｜为什么反演会有很多答案？](/posts/2026/10/mt-course-07a-inverse/)：零空间、存在性、唯一性与稳定性。
2. [02｜从残差到加权 RMS](/posts/2026/10/mt-course-07b-misfit/)：误差和数据失配。
3. [03｜模型粗糙度与正则化](/posts/2026/10/mt-course-07c-regularization/)：Occam 的约束思想。
4. [04｜雅可比矩阵到底是什么？](/posts/2026/10/mt-course-07d-jacobian/)：从数字表到模型敏感性。
5. [05｜梯度、Hessian 与 Gauss–Newton](/posts/2026/10/mt-course-07e-gauss-newton/)：如何计算模型更新。
6. [06｜L 曲线与恢复试验](/posts/2026/10/mt-course-07f-diagnostics/)：如何判断模型是否可信。

建议从第一篇开始；下方保留原先的浓缩概览，便于学完后整体复习。

---


## 1. 一个不能忽略的事实：反演可能有许多解

![正演与反演](/images/posts/mt-course/05-forward-inverse.svg)

**正演**：给定地下模型 $\mathbf m$（例如每层电阻率与厚度），使用物理规律计算预测观测 $\mathbf d_{\rm pred}=F(\mathbf m)$。**反演**：已知观测 $\mathbf d_{\rm obs}$，寻找使预测与观测足够接近的模型。

为什么“反过来求”突然变难？地表只能观测有限频率、有限站点的数据，有噪声，而地下有很多未知数。不同结构可以产生近似相同的地表响应，例如较薄且更导电的层与较厚但导电性稍弱的层，可能具有相近的积分效应。问题既存在由于数据有限导致的不足，也存在某些物理问题固有的非唯一性。

数学上，如果存在 $\Delta\mathbf m\ne0$ 但 $\mathbf G\Delta\mathbf m=0$，那么对于线性模型 $\mathbf d=\mathbf G\mathbf m$，$\mathbf m$ 与 $\mathbf m+\Delta\mathbf m$ 预测完全相同。这类无法被数据区分的方向属于**零空间**。零空间是反演中必须面对的结构，而不是计算机不够快造成的。

## 2. 先学最简单的一维线性拟合

假设观测满足 $y_j=ax_j+b+\varepsilon_j$。给定一些 $(x_j,y_j)$，寻找直线的斜率 $a$ 与截距 $b$。通常希望让

$$
\sum_{j=1}^N[y_j-(ax_j+b)]^2
$$

尽可能小。这就是最小二乘原理：测量值与预测值的差叫**残差**，残差平方和小表示拟合较好。MT 把这里的 $(a,b)$ 替换为成百上千个地下电阻率单元，把预测函数换成麦克斯韦正演 $F(\mathbf m)$。

## 3. 观测有误差，为什么要加权？

若第 $j$ 项数据的不确定度标准差是 $s_j$，应关注残差占“自身噪声水平”的比例：

$$
r_j(\mathbf m)=\frac{F_j(\mathbf m)-d_j}{s_j}.
$$

归一化数据失配函数常写作

$$
\Phi_d(\mathbf m)=\sum_{j=1}^N\left|\frac{F_j(\mathbf m)-d_j}{s_j}\right|^2.
$$

若处理阻抗的实部、虚部或视电阻率、相位，必须考虑它们的单位、误差定义与数据协方差。更通用的矩阵形式：

$$
\Phi_d=(F(\mathbf m)-\mathbf d)^T\mathbf C_d^{-1}(F(\mathbf m)-\mathbf d),
$$

这里为实数展开形式；若直接使用复数数据，需采用相应的 Hermitian（二次型共轭转置）写法。$\mathbf C_d$ 是数据误差协方差矩阵。误差估计不可靠时，即使“拟合到 RMS=1”也未必意味着模型正确。

## 4. 为什么还要“惩罚”模型？

如果单纯追求 $\Phi_d$ 最小，参数多时模型可能出现大量剧烈振荡，试图拟合观测噪声。于是增加一个反映模型平滑、简洁或偏离先验程度的约束项：

$$
\boxed{\Phi(\mathbf m)=\Phi_d(\mathbf m)+\lambda\Phi_m(\mathbf m)}.
$$

$\lambda>0$ 是正则化权重；较大 $\lambda$ 通常更强调模型约束，较小 $\lambda$ 更强调数据拟合。比如

$$
\Phi_m=\|\mathbf L(\mathbf m-\mathbf m_{\rm ref})\|_2^2.
$$

$\mathbf L$ 可为相邻单元差分算子，例如一维模型

$$
\mathbf L\mathbf m=
\begin{pmatrix}
m_2-m_1\\m_3-m_2\\\vdots\\m_n-m_{n-1}
\end{pmatrix}.
$$

这很容易理解：相邻两层参数差得越大，惩罚通常越大。常用 $m_j=\log\rho_j$ 或 $\log\sigma_j$ 使模型参数便于保持正值并允许跨数量级变化，具体取决于程序设定。

## 5. 反演为什么需要雅可比矩阵？

把“略微改变模型参数，预测数据会改变多少”写成偏导数。定义

$$
\boxed{J_{jk}=\frac{\partial F_j(\mathbf m)}{\partial m_k}}.
$$

$\mathbf J$ 称为**雅可比（灵敏度）矩阵**。第 $j$ 行对应某一条数据，第 $k$ 列对应地下某个模型参数。例如若调高第二层电阻率 1%，观测的某个频率相位变化多少，就反映在相应灵敏度里。实际问题可以用小扰动有限差分、伴随法等求灵敏度。

我们把非线性函数在当前模型 $\mathbf m$ 附近做一阶泰勒近似：

$$
F(\mathbf m+\Delta\mathbf m)\approx F(\mathbf m)+\mathbf J\Delta\mathbf m.
$$

这句话是非线性反演数值求解的基础：虽然完整问题不是线性的，但在足够小的邻域里可以近似为线性问题，求一次更新，再重新计算新的预测数据。

## 6. 梯度、Hessian 与 Newton 方法

对标量目标函数 $\Phi(\mathbf m)$，梯度 $\mathbf g=\nabla\Phi$ 是各参数偏导数组成的向量，指出局部上升最快的方向；Hessian $\mathbf H=\nabla\nabla\Phi$ 是二阶偏导数组成的矩阵，描述局部曲率。泰勒展开到二阶：

$$
\Phi(\mathbf m+\Delta\mathbf m)\approx
\Phi(\mathbf m)+\mathbf g^T\Delta\mathbf m+
\frac12\Delta\mathbf m^T\mathbf H\Delta\mathbf m.
$$

令近似式对 $\Delta\mathbf m$ 的导数为零：

$$
\mathbf H\Delta\mathbf m=-\mathbf g.
$$

这就是牛顿型方法的核心更新方程。实际地球物理反演可能存在病态矩阵，因此需要正则化、预条件或迭代线性求解器，不能粗暴地显式求逆。

## 7. Gauss–Newton 更新如何推导？

为了便于推导，先取

$$
\Phi=\frac12\|\mathbf W_d(F(\mathbf m)-\mathbf d)\|^2+
\frac{\lambda}{2}\|\mathbf L(\mathbf m-\mathbf m_{\rm ref})\|^2.
$$

令残差 $\mathbf r=F(\mathbf m)-\mathbf d$，对 $F$ 作一阶近似，再令增量二次函数最小，可得到

$$
\boxed{
(\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf J
+\lambda\mathbf L^T\mathbf L)\Delta\mathbf m
=
-\mathbf J^T\mathbf W_d^T\mathbf W_d\mathbf r
-\lambda\mathbf L^T\mathbf L(\mathbf m-\mathbf m_{\rm ref})
}.
$$

左边括号内的矩阵可看作正则化后的近似 Hessian。相比完整 Newton 方法，Gauss–Newton 忽略了由正演函数二阶导数产生的一部分项，计算通常更可行。

完成一次更新后，不一定直接走完整步：

$$
\mathbf m_{\rm next}=\mathbf m+\gamma\Delta\mathbf m.
$$

$\gamma$ 是步长，通常借助**线搜索**选择，以避免一步走得太远，反而使目标函数上升。

## 8. Occam 反演为什么强调“尽量简单”？

当许多地下模型都能解释数据时，Occam 思想选择在满足目标数据拟合水平的候选模型中尽量平滑或粗糙度较小的模型。通过约束模型粗糙度，并选择适当的目标误差，可以得到没有过多不必要细节的结构。

它的优点是减少人为加入尖锐界面的风险；局限是**真实地下若确实存在断层、薄导电层或突变，平滑偏好也可能把它们抹宽或弱化**。因此 Occam 模型是一种受数据与先验共同约束的结果，不能把每个颜色边界理解成真实地质界面。

本章还介绍 L 曲线等权衡工具：把模型粗糙度与数据拟合误差画在一起，寻找折中。需要强调：曲线上的一个“拐角”也并不是绝对真实模型的证明。

## 9. 实际结果需要做什么验证？

- **数据拟合**：逐频率、逐站点比较观测与预测的视电阻率、相位、阻抗，避免只看一个全局 RMS。
- **起始模型与先验**：尝试不同起始模型或先验结构，检查较深结构是否稳定。
- **灵敏度与分辨率**：对最关心的异常作扰动或恢复实验，不要把平滑彩图等同于同等分辨率。
- **异常成因**：即使低阻结构稳定，仍需岩石物性、温压、流体、地震以及构造背景来判断原因。
- **模型复杂度**：比较一维、二维、三维假设是否与数据维度分析一致。

具体软件可因数据格式和算法有所差别，但遵守相同的正演—失配—更新—检验流程。

### 小练习

1. $\Phi=(m-3)^2$ 的梯度是什么？答：$2(m-3)$；最小值在 $m=3$。
2. 如果增大 $\lambda$，通常更强调什么？答：模型约束项，常产生更平滑或更接近先验的结果。
3. “模拟数据可以和观测完全重合”是否等于“地下结构唯一且真实”？答：不等于，仍可能有零空间、系统误差和先验偏好 **下一讲：**[从电阻率结构到地质解释](/posts/2026/10/mt-course-08-interpretation/)。