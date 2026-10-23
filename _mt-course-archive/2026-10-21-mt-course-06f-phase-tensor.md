---
title: "数据分析 06｜相位张量怎么来的：先用数字理解矩阵求逆"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-06f-phase-tensor/
excerpt: "由 Z=X+iY 开始，手算二阶矩阵的逆及乘法，再完整推导 Φ=X⁻¹Y 的畸变不变性和椭圆表示。"
tags: [大地电磁, 相位张量, 矩阵求逆, 数学基础]
---
{% include mt-course-math.html %}

我们在上篇见到一个麻烦：观测到的阻抗可能是

$$
\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg},
$$

其中 $\mathbf C$ 表示近地表的实数电流性畸变矩阵。能不能构造某个数学量，利用矩阵运算让 $\mathbf C$ 自动抵消？

这正是**相位张量（phase tensor）**的重要动机。为了真正理解它，我们先从最容易算的二阶矩阵开始。

## 一、为什么要把复数矩阵分开？

四个复阻抗元素都可以写成“实部+虚部”。所以整个矩阵可记作

$$
\boxed{\mathbf Z=\mathbf X+i\mathbf Y}.
$$

$\mathbf X$ 表示四个元素实部组成的实矩阵，$\mathbf Y$ 表示虚部组成的实矩阵。这里的 $X,Y$ 是**两个完整矩阵**，不是横坐标 $x,y$；符号相近但意义不同。

例如

$$
\mathbf Z=
\begin{pmatrix}1+2i&3+i\\2-i&4+3i\end{pmatrix},
$$

则

$$
\mathbf X=\begin{pmatrix}1&3\\2&4\end{pmatrix},
\quad
\mathbf Y=\begin{pmatrix}2&1\\-1&3\end{pmatrix}.
$$

## 二、矩阵的“倒数”究竟是什么？

普通数字里 $3\times(1/3)=1$。矩阵中对应的概念是**逆矩阵**：

$$
\mathbf X^{-1}\mathbf X=\mathbf I,
$$

其中 $\mathbf I=\begin{pmatrix}1&0\\0&1\end{pmatrix}$。不是什么矩阵都有逆；例如两行完全成比例，逆矩阵便不存在。

对

$$
\mathbf X=\begin{pmatrix}a&b\\c&d\end{pmatrix},
$$

若行列式 $ad-bc\ne0$，逆矩阵为

$$
\boxed{\mathbf X^{-1}=\frac1{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}.}
$$

这条公式可以靠乘法验证。比如

$$
\mathbf X=\begin{pmatrix}2&0\\0&4\end{pmatrix},
\quad
\mathbf X^{-1}=\begin{pmatrix}1/2&0\\0&1/4\end{pmatrix}.
$$

相乘确实是单位矩阵。它就像分别把两个方向上放大的 2 倍与 4 倍缩回去。

## 三、亲手算一个相位张量

定义

$$
\boxed{\boldsymbol\Phi=\mathbf X^{-1}\mathbf Y}.
$$

假设阻抗实部、虚部为

$$
\mathbf X=\begin{pmatrix}2&0\\0&4\end{pmatrix},
\quad
\mathbf Y=\begin{pmatrix}1&0\\0&2\end{pmatrix}.
$$

因为 $\mathbf X^{-1}=\begin{pmatrix}1/2&0\\0&1/4\end{pmatrix}$，所以

$$
\boldsymbol\Phi=
\begin{pmatrix}1/2&0\\0&1/4\end{pmatrix}
\begin{pmatrix}1&0\\0&2\end{pmatrix}
=
\begin{pmatrix}1/2&0\\0&1/2\end{pmatrix}.
$$

它是一个实数矩阵，并且四个阻抗复数中的幅值信息已经被一定程度上消去。

## 四、最重要的证明：为什么畸变会抵消？

我们假设只有近地表实数电流性电场畸变，且磁场不受对应畸变：

$$
\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg}.
$$

把区域阻抗写为

$$
\mathbf Z_{\rm reg}=\mathbf X+i\mathbf Y.
$$

因为 $\mathbf C$ 是实矩阵，所以

$$
\mathbf Z_{\rm obs}
=\mathbf C\mathbf X+i\mathbf C\mathbf Y.
$$

这说明观测实部矩阵、虚部矩阵分别是

$$
\mathbf X_{\rm obs}=\mathbf C\mathbf X,
\quad
\mathbf Y_{\rm obs}=\mathbf C\mathbf Y.
$$

按相位张量定义：

$$
\boldsymbol\Phi_{\rm obs}
=(\mathbf C\mathbf X)^{-1}(\mathbf C\mathbf Y).
$$

我们需要一个关于矩阵求逆的规则：

$$
(\mathbf A\mathbf B)^{-1}=\mathbf B^{-1}\mathbf A^{-1}.
$$

为什么顺序要反过来？因为先做 A 再做 B 的整体操作，要恢复它时必须逆着运算顺序撤销；可以直接将右式乘回去验证。

因此

$$
\begin{aligned}
\boldsymbol\Phi_{\rm obs}
&=\mathbf X^{-1}\mathbf C^{-1}\mathbf C\mathbf Y\\
&=\mathbf X^{-1}\mathbf I\mathbf Y\\
&=\boxed{\boldsymbol\Phi_{\rm reg}}.
\end{aligned}
$$

看到了吗？中间的 $\mathbf C^{-1}\mathbf C$ 等于单位矩阵，于是某一类电流性畸变被消去了。

这个结论的条件非常重要：$\mathbf C$ 与 $\mathbf X$ 要可逆，$\mathbf C$ 要在研究频带内适合用实矩阵描述，磁场畸变可忽略，并且其他更复杂的频率相关感应影响不能直接等同于这一模型。

## 五、相位张量元素怎么实际计算？

令

$$
\mathbf X=\begin{pmatrix}X_{xx}&X_{xy}\\X_{yx}&X_{yy}\end{pmatrix},
\quad
\mathbf Y=\begin{pmatrix}Y_{xx}&Y_{xy}\\Y_{yx}&Y_{yy}\end{pmatrix}.
$$

行列式是

$$
D=X_{xx}X_{yy}-X_{xy}X_{yx}.
$$

根据二阶逆矩阵公式：

$$
\boldsymbol\Phi=\frac1D
\begin{pmatrix}
X_{yy}Y_{xx}-X_{xy}Y_{yx}&X_{yy}Y_{xy}-X_{xy}Y_{yy}\\
X_{xx}Y_{yx}-X_{yx}Y_{xx}&X_{xx}Y_{yy}-X_{yx}Y_{xy}
\end{pmatrix}.
$$

当 $D$ 接近零时，求逆会非常不稳定，噪声可被显著放大。这是实际相位张量估计必须检查的数值问题。

## 六、为什么要把相位张量画成椭圆？

![相位张量的圆和椭圆示意](/images/posts/mt-course/25-phase-tensor-ellipse.svg)

一个二维矩阵可以描述不同方向上的变换。对相位张量常采用主轴、椭圆率和偏斜角绘制椭圆：椭圆长短轴反映不同主方向的响应差异，方向角描述其方位，颜色还可表示偏斜或其他指标。

对于理想一维各向同性结构，相位张量与单位矩阵成比例，因此通常画成圆。理想二维模型在走向坐标中没有偏斜，主方向可以与地质走向或其垂直方向联系起来。但**画出椭圆不等于证明地下必然二维，更不保证它的长轴就是唯一正确的地质走向**。

## 七、四个值得认识的数

记

$$
\boldsymbol\Phi=
\begin{pmatrix}a&b\\c&d\end{pmatrix}.
$$

可以构造

$$
\Phi_1=\frac{a+d}{2},\quad
\Phi_2=ad-bc,\quad
\Phi_3=\frac{b-c}{2}.
$$

它们分别涉及迹、行列式和反对称部分。常用的方向参数为

$$
\alpha=\frac12\operatorname{atan2}(b+c,a-d),
$$

$$
\beta=\frac12\operatorname{atan2}(b-c,a+d).
$$

$\operatorname{atan2}(y,x)$ 与普通 $\arctan(y/x)$ 不同，它能利用分子与分母的符号处理象限。$\beta$ 常叫相位张量偏斜角，理想二维模型有 $\beta=0$，一般三维结构可能得到非零值；但噪声及接近退化的矩阵会使角度估计不稳定。

还可定义归一化各向异性指标

$$
\lambda=
\frac{\sqrt{(a-d)^2+(b+c)^2}}
{\sqrt{(a+d)^2+(b-c)^2}},
$$

分母须非零。它与 $\beta$ 可共同辅助维度判断：理想一维两者都为零，理想二维一般 $\beta=0$ 但 $\lambda$ 可以非零；不能把这些理想条件无条件套在有噪声数据上。

## 八、相位张量不能做什么？

它只有四个实数分量，而原始复阻抗有八个实数分量。因此**相位张量不包含阻抗的完整振幅信息**。使用它诊断电性维数和走向很有用，但不能认为只凭这四个数就能恢复原来的完整阻抗张量，也不能直接替代完整幅相反演。

**练习：** 若 $\mathbf X=2\mathbf I$、$\mathbf Y=3\mathbf I$，则 $\boldsymbol\Phi=(3/2)\mathbf I$，对应理想的各向同性相位张量形态。

[下一节：走向估计与维度诊断](/posts/2026/10/mt-course-06g-strike/) · [专题目录](/posts/2026/10/mt-course-06-analysis/)