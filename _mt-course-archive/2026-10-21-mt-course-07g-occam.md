---
title: "地球物理反演 07｜Occam 为什么偏爱平滑模型：从相邻两层之差开始"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-07g-occam/
excerpt: "手算模型粗糙度，理解正则化、参考模型、目标 RMS、L 曲线，以及 Occam 线性化更新与 Gauss–Newton 的联系。"
tags: [大地电磁, Occam反演, 正则化, 平滑约束]
---
{% include mt-course-math.html %}

我们终于知道怎样用 Gauss–Newton 更新地下电阻率模型。但仍有一个根本困难没有解决：**如果好几种地下模型都能解释观测，计算机究竟应该选哪一种？**

Occam（奥卡姆）反演给出了一个很有影响力的选择原则：在达到合理数据拟合水平的模型中，选择相对平滑、没有不必要起伏的结果。

这听起来像一句哲学口号，其实可以落实到非常具体的数学运算。

## 一、怎样量化模型到底有多“粗糙”？

先只考虑三层地下，模型参数分别为 $m_1,m_2,m_3$。注意实际反演里常选择 $m_j=\log\rho_j$ 或 $\log\sigma_j$，而不是直接对电阻率本身求差。

第一种模型：

$$\mathbf m_A=\begin{pmatrix}1\\3\\1\end{pmatrix}.$$

相邻层的参数差为

$$m_2-m_1=2,\qquad m_3-m_2=-2.$$

把它们分别平方，再求和：

$$
\Phi_m(\mathbf m_A)=2^2+(-2)^2=8.
$$

第二种模型：

$$\mathbf m_B=\begin{pmatrix}1\\2\\3\end{pmatrix}.$$

相邻参数差为 1 和 1，所以

$$
\Phi_m(\mathbf m_B)=1^2+1^2=2.
$$

第二种模型的参数变化比较平缓，因此我们定义的粗糙度更小。

![模型粗糙度示意](/images/posts/mt-course/32-roughness.svg)

这里并没有判断哪一个模型更真实；我们只是定义了一个可以计算的“平滑程度”。

## 二、怎样用一个矩阵求所有相邻差？

定义差分矩阵

$$
\mathbf L=
\begin{pmatrix}
-1&1&0\\
0&-1&1
\end{pmatrix}.
$$

将它乘以三层模型：

$$
\mathbf L\mathbf m=
\begin{pmatrix}
m_2-m_1\\
m_3-m_2
\end{pmatrix}.
$$

因此粗糙度函数就可以写成

$$
\boxed{\Phi_m(\mathbf m)=\|\mathbf L\mathbf m\|^2}.
$$

$\|\cdot\|^2$ 表示所有元素平方再求和。更一般地，如果希望模型接近某个参考结构 $\mathbf m_{\rm ref}$，也可以写

$$
\boxed{
\Phi_m=\|\mathbf L(\mathbf m-\mathbf m_{\rm ref})\|^2.
}
$$

对于非均匀厚度网格，差分还应考虑层厚和深度间距。二维、三维模型也可以对水平与垂直方向采用不同的平滑权重。

## 三、为什么还需要数据失配？

如果只让粗糙度最小，一份所有层参数完全相等的常数模型就可能得到 $\Phi_m=0$。但它未必能解释真实测量。

所以我们既要数据拟合，也要模型平滑。常见的正则化目标写为

$$
\boxed{
\Phi(\mathbf m)=\Phi_d(\mathbf m)+\lambda\Phi_m(\mathbf m).
}
$$

这里 $\Phi_d$ 是误差标准化后的数据失配，$\Phi_m$ 是粗糙度，$\lambda>0$ 是正则化权重。

如果 $\lambda$ 很大，程序更看重平滑，可能把真实异常都抹掉；如果 $\lambda$ 很小，程序更看重拟合，可能产生很多不必要的起伏。

**请留意不同软件、文章采用的权重约定。** 有些写成 $\Phi_m+\lambda^{-1}\Phi_d$，此时 $\lambda$ 和上式的直观大小关系需要重新解释，不能只凭一个希腊字母比较参数大小。

## 四、Occam 的目标与普通加权相加有什么不同？

Occam 的经典思路常表述为一个带约束的优化问题：

$$
\boxed{
\text{使 }\Phi_m(\mathbf m)\text{ 最小，条件是 }
\Phi_d(\mathbf m)\leq\Phi_d^*.
}
$$

其中 $\Phi_d^*$ 是事先设定的合理目标失配。也就是说：先规定“预测与观测要接近到什么程度”，然后在满足这一条件的模型中找尽可能平滑的那个。

在数值求解时，可以采用拉格朗日乘子或者加权目标函数

$$
\Phi_m+\mu\Phi_d
$$

把问题转化为便于优化的形式，调节 $\mu$ 以寻找满足目标拟合的解。这里的 $\mu$ 不应与前面电磁学磁导率 $\mu_0$ 混淆。

这也解释了为什么 Occam 和 Gauss–Newton 不是“只能二选一”的两种不相干技术：**Occam 决定优化偏好与目标拟合原则，具体迭代仍可以用线性化和类似 Gauss–Newton 的数值计算来完成。**

## 五、Occam 怎么得到下一次模型？

设当前模型为 $\mathbf m_k$，希望得到下一次模型 $\mathbf m_{k+1}$。用雅可比做局部线性化：

$$
F(\mathbf m_{k+1})
\approx F(\mathbf m_k)+
\mathbf J_k(\mathbf m_{k+1}-\mathbf m_k).
$$

为了整理公式，定义

$$
\widehat{\mathbf d}_k=
\mathbf d-F(\mathbf m_k)+\mathbf J_k\mathbf m_k.
$$

于是当前线性近似下的观测方程可以写成

$$
\widehat{\mathbf d}_k\approx\mathbf J_k\mathbf m_{k+1}.
$$

假设模型约束采用正定的 $\mathbf C_m^{-1}$，数据权重采用 $\mathbf C_d^{-1}$，构造

$$
\Phi_{\rm approx}(\mathbf m_{k+1})
=
(\mathbf m_{k+1}-\mathbf m_{\rm ref})^T
\mathbf C_m^{-1}
(\mathbf m_{k+1}-\mathbf m_{\rm ref})
+\mu(\widehat{\mathbf d}_k-\mathbf J_k\mathbf m_{k+1})^T
\mathbf C_d^{-1}
(\widehat{\mathbf d}_k-\mathbf J_k\mathbf m_{k+1}).
$$

第一项管模型形态，第二项管模型预测与数据的距离。将目标对 $\mathbf m_{k+1}$ 求导并令其为零，整理得到

$$
\boxed{
\begin{aligned}
(&\mathbf C_m^{-1}
+\mu\mathbf J_k^T\mathbf C_d^{-1}\mathbf J_k)
\mathbf m_{k+1}
={}&\mathbf C_m^{-1}\mathbf m_{\rm ref}\\
&+\mu\mathbf J_k^T\mathbf C_d^{-1}\widehat{\mathbf d}_k.
\end{aligned}
}
$$

如果左边矩阵可逆，可以形式上写成求逆后的解。实际计算依然更倾向于使用稳定的线性方程求解算法。

这段推导为什么要放在这里？是为了看见它与 Gauss–Newton 的共性：**两者都使用当前模型的雅可比矩阵，把非线性反演暂时近似为一个线性系统**。主要区别体现在目标函数、约束形式及权重选择原则上。

严格的“相邻差分平滑”对应 $\mathbf L^T\mathbf L$，它可能是半正定而非可逆的；使用上面的 $\mathbf C_m^{-1}$ 写法时，需要有适当正定先验或额外约束来保证可解，不能机械地把一个奇异差分算子求逆。

## 六、为什么 RMS 越低，模型反而可能越不可信？

![拟合与复杂度之间的 L 曲线示意](/images/posts/mt-course/33-lcurve.svg)

假设模型 A 的 nRMS 是 1.1，模型 B 的 nRMS 是 1.0。两者都处于合理拟合范围，但模型 B 为达到这一点增加了一大堆高度起伏的低阻、高阻小块。

这些新出现的复杂结构可能只是在拟合观测噪声。数据拟合从 1.1 改善到 1.0，并不足以单独证明新增异常真实存在。

有时可以把不同权重得到的**数据失配**与**模型粗糙度**画在一张图上，形成所谓的 **L 曲线（L-curve）**。拐角区域常被用于寻找折中，但它也会因噪声、尺度、模型参数化而变化，不能把拐角当成真实地质结构的数学证明。

## 七、Occam 的平滑结果是否一定接近真实地下？

不一定。真实地下完全可能有尖锐断层、突然改变的地层界面、很薄的导电层。平滑约束可能将这些结构的边界变宽，甚至把多个相邻异常合并。

因此 Occam 应被理解为：**当数据不足以证明某些尖锐细节时，先选择比较简单、平滑的解释。** 这个模型便于看清数据稳定要求的主要异常，但仍需进行灵敏度、分辨率与恢复实验。

## 八、自检

**①** 两个模型同样能满足目标 RMS，为什么不能只选结构变化最大的那个？因为复杂起伏可能是数据并不要求的额外自由度。

**②** 为什么不能一味把正则化权重调得越来越大？因为模型可能过分平滑，无法再解释观测。

**③** 平滑反演中出现一个 10 km 宽的低阻区，是否说明真实导电地层一定宽 10 km？不能，范围受数据分辨率与平滑约束共同影响。

[下一节：反演结果怎样验证？](/posts/2026/10/mt-course-07h-validation/) · [反演专题目录](/posts/2026/10/mt-course-07-inversion/)