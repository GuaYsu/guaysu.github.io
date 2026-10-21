---
title: "MT 讲义 06｜阻抗张量、维度、相位张量与静态位移"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-06-analysis/
excerpt: "从二维矩阵坐标旋转和阻抗张量讲到 1D/2D/3D 识别、旋转不变量、Groom–Bailey 畸变、相位张量、走向与静态位移。"
tags: [大地电磁, 阻抗张量, 相位张量, 畸变分析]
---
{% include mt-course-math.html %}

## 1. 为什么不能一直用 $Z=E_x/H_y$？

前面讨论的是电性只随深度变化的理想一维模型。在野外，地下导电断层、盆地边界或三维异常可以使不同水平分量互相耦合，所以需要使用完整的**阻抗张量**：

$$
\begin{pmatrix}E_x\\E_y\end{pmatrix}=
\underbrace{\begin{pmatrix}
Z_{xx}&Z_{xy}\\Z_{yx}&Z_{yy}
\end{pmatrix}}_{\mathbf Z}
\begin{pmatrix}H_x\\H_y\end{pmatrix}.
$$

把矩阵乘法逐行写开：

$$
E_x=Z_{xx}H_x+Z_{xy}H_y,\qquad
E_y=Z_{yx}H_x+Z_{yy}H_y.
$$

每个 $Z_{ij}$ 一般都是随频率变化的复数。$Z_{xy}$ 指 $H_y$ 对 $E_x$ 的线性关系系数，**不是**“电场在 $x$ 位置与磁场在 $y$ 位置的比值”。若磁场用磁感应强度 $\mathbf B$ 而非 $\mathbf H$ 表示，响应系数还涉及 $\mathbf B=\mu_0\mathbf H$ 的单位换算；课件将对应的 MT 张量记作 $\mathbf M$，有 $\mathbf Z=\mu_0\mathbf M$。垂直磁场也可用 tipper 表示：

$$
H_z=T_{zx}H_x+T_{zy}H_y.
$$

它提供另一类空间维度与异常方向的信息，不能把 tipper 当成电阻率本身。

## 2. 1D、2D、3D 是怎样区分的？

先约定：这里的维度指**地下电性模型依赖几个空间坐标**。

- **1D**：$\rho=\rho(z)$，水平方向均匀；各向同性模型在任何水平旋转下有 $Z_{xx}=Z_{yy}=0$、$Z_{xy}=-Z_{yx}$，理想的垂直传递函数为 0。
- **2D**：例如 $\rho=\rho(y,z)$，沿走向 $x$ 不改变。坐标轴与走向一致时，阻抗的两个对角元为 0，但通常 $Z_{xy}\ne-Z_{yx}$；前者后者对应不同极化模式。
- **3D**：$\rho=\rho(x,y,z)$，一般四个张量元素均可能非零，tipper 也可能两方向均非零。

这里强调“理想模型”很重要：实际观测有噪声，仪器方向有误差，浅部畸变会让简洁的零元素不再等于严格的零。仅凭某一个零值或某个频点，不能宣布地下恰好一维或二维。

## 3. TE 和 TM：与走向有关的两种模式

对于沿 $x$ 走向不变的二维结构：

**TE（横向电）模式**：电场平行走向，电流沿走向流动，通常考虑 $(E_x,H_y)$；**TM（横向磁）模式**：磁场平行走向，电流横穿走向，通常考虑 $(E_y,H_x)$。不同模式对断层、盆地边界等二维结构的敏感性不同。

一个很常见的错误是把“仪器 $x$ 轴”自动当成“地质走向”。如果仪器布设方向与真正地质走向不一致，必须先旋转张量，不能直接拿未经旋转的 $Z_{xy}$ 当 TE 响应。

## 4. 为什么旋转坐标会出现两个矩阵？

假设把测量坐标转过角度 $\theta$，按课件约定定义

$$
\mathbf R=
\begin{pmatrix}
\cos\theta&\sin\theta\\
-\sin\theta&\cos\theta
\end{pmatrix}.
$$

新坐标里的电场、磁场分别为 $\mathbf E'=\mathbf R\mathbf E$、$\mathbf H'=\mathbf R\mathbf H$。原坐标下 $\mathbf E=\mathbf Z\mathbf H$，因此

$$
\mathbf E'=\mathbf R\mathbf Z\mathbf H
=\mathbf R\mathbf Z\mathbf R^T\mathbf H'.
$$

与 $\mathbf E'=\mathbf Z'\mathbf H'$ 比较：

$$
\boxed{\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T}.
$$

这里 $T$ 表示**转置**，因为旋转矩阵满足 $\mathbf R^{-1}=\mathbf R^T$。一个矩阵在左边是为了旋转输出电场，在右边是为了把输入磁场换回相同坐标。这个推导比死记公式更可靠。

### 简单例子

若旋转 $90^\circ$，$\cos90^\circ=0,\sin90^\circ=1$，则

$$
\mathbf R=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
$$

把它乘入上式，可以得到新坐标下的所有四个阻抗分量。无论如何旋转，同一个物理观测并没有改变，只是数字的表达方式变化了。

## 5. 旋转不变量有什么用？

如果某种组合不受坐标旋转影响，它可以帮助我们比较不同仪器方向的结果。比如迹、反迹与行列式：

$$
\operatorname{tr}\mathbf Z=Z_{xx}+Z_{yy},
$$

$$
\operatorname{atr}\mathbf Z=Z_{xy}-Z_{yx},
$$

$$
\det\mathbf Z=Z_{xx}Z_{yy}-Z_{xy}Z_{yx}.
$$

对二维正交旋转，迹与行列式保持不变，反迹在此固定旋转方向约定下也不变。需要留意，这几个量本身是否为零或满足特定关系，并不能独立解决所有畸变与维度判断问题。

课件进一步介绍 Weaver–Agarwal–Lilley（WAL）不变量 $I_1,\ldots,I_7$：其中一些描述不同水平极化的差别，另一些关联畸变与维度。它们的数学表达式较长，实际使用时要明确实部和虚部、归一化以及阈值如何定义，不能只背“1D 就都为零”。特别地，$I_7$ 在某些一维理想情况下可出现未定义情形。

## 6. 浅部小异常如何“扭曲”深部信息？

想象真正想研究的区域电场为 $\mathbf E_{\rm reg}$。近地表存在一个尺度很小的导电体，它可以改变局部电流的分布，使测到的水平电场发生混合：

$$
\mathbf E_{\rm obs}=\mathbf C\mathbf E_{\rm reg},
$$

其中 $\mathbf C$ 是描述近地表、近似实数和频率无关的**电流性（galvanic）畸变矩阵**。若假设磁场不发生相应畸变，则

$$
\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg}.
$$

这就是后续 Groom–Bailey（GB）畸变分解的重要出发点。实际畸变可能包含整体增益、扭转（twist）、剪切（shear）等分量，并与区域走向联立估计。单靠某一个频点、某一个站的八个实数阻抗观测，某些分量可能不唯一，所以课件讨论多个频率、多个站点的联合求解。

## 7. 相位张量：如何绕开某些畸变？

写 $\mathbf Z=\mathbf X+i\mathbf Y$，其中 $\mathbf X$、$\mathbf Y$ 分别是阻抗的实部矩阵、虚部矩阵。定义相位张量

$$
\boxed{\boldsymbol\Phi=\mathbf X^{-1}\mathbf Y}.
$$

为什么它对上述实数畸变不敏感？因为在假设 $\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg}$ 下，

$$
\mathbf X_{\rm obs}=\mathbf C\mathbf X_{\rm reg},\quad
\mathbf Y_{\rm obs}=\mathbf C\mathbf Y_{\rm reg},
$$

于是

$$
\begin{aligned}
\boldsymbol\Phi_{\rm obs}
&=(\mathbf C\mathbf X_{\rm reg})^{-1}(\mathbf C\mathbf Y_{\rm reg})\\
&=\mathbf X_{\rm reg}^{-1}\mathbf C^{-1}\mathbf C\mathbf Y_{\rm reg}\\
&=\boldsymbol\Phi_{\rm reg}.
\end{aligned}
$$

看到这里，矩阵逆为什么重要就清楚了：$\mathbf C^{-1}\mathbf C=\mathbf I$，它们抵消了。这个结论需要 $\mathbf C$ 和 $\mathbf X$ 可逆、仅有这种频率无关的电流性电场畸变、磁场不畸变等条件；其他畸变或噪声不保证消除。

相位张量可用椭圆表示；主轴方向、椭圆率和偏斜角 $\beta$ 常用于分析维度与走向。**相位张量舍弃了阻抗整体振幅信息**，因此不能仅凭它独立完成通常需要绝对电阻率尺度的反演。

## 8. 从维度判断到走向判断

课件列举经典 Swift skew、Bahr phase-sensitive skew、WAL 不变量、极坐标图（peanut diagrams）、相位张量椭圆、Swift/Bahr 走向估计等方法。Swift skew 曾用约 0.2 等经验阈值区分较简单的维度响应，但它对噪声及电流性畸变敏感，阈值只能作为筛查线索。

常见工作顺序是：分别绘制不同周期的视电阻率与相位、skew 与 tipper、相位张量椭圆，检查变化是否稳定；确定有意义的近二维频带之后，再估计走向并把张量旋转到相应坐标中。不同周期探测深度不同，所以“走向随周期改变”本身可能具有地质含义，也可能来自噪声或畸变，不能先入为主。

## 9. 静态位移（static shift）

近地表电流性畸变可使测点阻抗振幅近似被一个实数比例因子 $g$ 改变：

$$
Z_{\rm obs}(\omega)\approx gZ_{\rm reg}(\omega).
$$

因此视电阻率会整体平移：

$$
\rho_{a,\rm obs}\approx g^2\rho_{a,\rm reg},
$$

但当 $g$ 为正实数、频率无关时，相位几乎不变。这就是为什么两条视电阻率曲线可能形状相似却在对数纵轴上上下错开。课件讨论与附近测点、TEM 或其他地球物理数据联合约束静态位移；纯靠曲线外观任意“拉齐”有风险，因为真实电阻率差异也会造成变化。

课程还介绍用 MTpy 绘制阻抗、视电阻率与相位，检查相位张量与电磁响应。程序只是工具，**诊断前提是否成立**才是最重要的事。

### 自检

1. 两个地层模型都产生 $Z_{xx}\ne0$，是否都一定是三维？答：不一定，坐标未对准走向或畸变也会使对角元素非零。
2. 旋转张量为什么要 $\mathbf R\mathbf Z\mathbf R^T$？答：输入与输出矢量都必须同时变换。
3. 相位张量能否完全替代原始阻抗？答：不能，它丢失整体幅值信息。
4. 静态位移能否单独由相位曲线纠正？答：一般不能，往往需要额外约束。

**课件定位：** Ch.5 第 4–68 页，主要推导第 5–12、16–32 页，诊断方法与实际例子第 33–68 页。**下一讲：**[反演与 Occam](/posts/2026/10/mt-course-07-inversion/)。