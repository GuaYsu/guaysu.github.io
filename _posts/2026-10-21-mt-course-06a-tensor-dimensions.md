---
title: "09｜阻抗张量、地质维度与坐标旋转｜合并完整版"
date: 2026-10-21
modified: 2026-10-23
permalink: /posts/2026/10/mt-course-06a-tensor-dimensions/
excerpt: "保留每节数学推导、定义、实际分析方法与自检练习。"
tags: [大地电磁, 阻抗张量, 数学基础]
---
{% include mt-course-math.html %}

完整合并原有专题，目录可跳转到各原章节。

## 本篇目录
- [阻抗张量与数据分析｜从矩阵入门到相位张量、走向和静态位移](#sec-06-analysis)
- [数据分析 01｜阻抗张量为什么是四个复数：从矩阵乘法开始](#sec-06a-tensor)
- [数据分析 02｜一维、二维、三维怎么判断？TE 和 TM 又是什么？](#sec-06b-dimensionality)
- [数据分析 03｜为什么旋转阻抗张量，要左右各乘一个矩阵？](#sec-06c-rotation)


---

<a id="sec-06-analysis"></a>
## 第1单元｜阻抗张量与数据分析｜从矩阵入门到相位张量、走向和静态位移

上一专题，我们已经从五路时间序列提取傅里叶系数，并估计出每个频率的复数阻抗张量。现在真正有意思的问题出现了：**如何从四个复数判断地下是怎样的结构？测量方向是否正确？近地表的小型异常是否改变了我们想研究的区域响应？**

这一专题把曾经挤在一篇文章里的内容拆成八个逐步深入的小单元。先用数字和矩阵乘法做简单练习，再走向旋转不变量、相位张量与实际诊断。

![阻抗张量把两个磁场分量映射为两个电场分量](/images/posts/mt-course/21-impedance-matrix.svg)

### 八个小单元

1. **[阻抗张量为什么是四个复数？](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06a-tensor)**  
   从 $E_x=ZH_y$ 推广到两个输入、两个输出，手算矩阵乘法，认识 tipper 与 $B/H$ 单位差异。

2. **[一维、二维、三维怎么判断？TE 和 TM 又是什么？](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06b-dimensionality)**  
   地下维数的真正含义、二维构造走向、两种极化，以及为什么对角元不为零也未必意味着三维。

3. **[为什么张量旋转要左右各乘一个矩阵？](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06c-rotation)**  
   从向量投影、三角函数与转置入手，完整推导 $\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T$，再介绍迹、反迹与行列式。

4. **[WAL 不变量与维度判断](/posts/2026/10/mt-course-06d-invariants/)**  
   先把四个复阻抗拆成八个实数，再组合旋转不变量，认识 $I_1$ 至 $I_7$、$Q$ 和经典 Swift skew。

5. **[Groom–Bailey 畸变从哪里来？](/posts/2026/10/mt-course-06e-distortion/)**  
   浅部导体、界面电荷、twist、shear、anisotropy 与 site gain 的物理意义及不可辨识问题。

6. **[相位张量：为什么求逆后畸变能抵消？](/posts/2026/10/mt-course-06f-phase-tensor/)**  
   手算二阶逆矩阵，推导 $\boldsymbol\Phi=\mathbf X^{-1}\mathbf Y$，理解相位张量椭圆及偏斜角。

7. **[如何估计地质走向？](/posts/2026/10/mt-course-06g-strike/)**  
   Swift/Bahr/WAL、peanut 极坐标图与相位张量方位，认识角度歧义、噪声与频率变化。

8. **[静态位移为什么会骗过反演？](/posts/2026/10/mt-course-06h-static-shift/)**  
   由 $Z_{\rm obs}=gZ_{\rm reg}$ 推导 $\rho_{a,\rm obs}=g^2\rho_{a,\rm reg}$，比较 TE/TM 畸变，学习诊断与 MTpy 工作流。

### 为什么这些内容要在反演之前学？

反演程序只知道我们交给它什么阻抗、相位和误差。它不会自动知道某个测点的电极方向写错了，或某些视电阻率因浅部异常整体抬高了四倍。若输入数据的解释前提有问题，反演结果可以拟合得非常漂亮，却形成不合理的地下电阻率结构。

所以张量分析的任务是：先明确数据代表什么、与怎样的地下维度假设相容，以及哪些干扰和系统性畸变需要纳入不确定度或额外约束。

到这里学会的复数、矩阵、旋转、偏导与统计误差思想，在[下一专题：反演与 Occam](/posts/2026/10/mt-course-07-inversion/)都会继续发挥作用。

[返回全系列目录](/posts/2026/10/mt-course/)。


---

<a id="sec-06a-tensor"></a>
## 第2单元｜数据分析 01｜阻抗张量为什么是四个复数：从矩阵乘法开始

在一维水平层状模型中，一个 $Z_{xy}=E_x/H_y$ 就能描述一种极化的响应。地下具有横向变化后，每个水平电场分量都可能同时受到两个水平磁场分量的影响，于是需要四个复数系数。

### 一、先从一个输入和一个输出开始

设一个频率上的复数磁场是 $H_y$，电场是 $E_x$，理想模型满足

$$E_x=Z_{xy}H_y.$$

这只是“输出=比例系数×输入”。例如 $H_y=2$、$Z_{xy}=3+4i$，则 $E_x=6+8i$。复数系数同时包含振幅比例和相位关系。

现在地下可能沿水平方向变化。$E_x$ 不一定只由 $H_y$ 控制，也可能与 $H_x$ 有关。因此用两个系数：

$$E_x=Z_{xx}H_x+Z_{xy}H_y.$$

同理，另一个输出 $E_y$ 也需要两个系数：

$$E_y=Z_{yx}H_x+Z_{yy}H_y.$$

两行叠在一起，就成了阻抗张量。

![从两路磁场到两路电场](/images/posts/mt-course/21-impedance-matrix.svg)

### 二、矩阵只是两行方程的简写

$$
\boxed{
\begin{pmatrix}E_x\\E_y\end{pmatrix}=
\begin{pmatrix}Z_{xx}&Z_{xy}\\Z_{yx}&Z_{yy}\end{pmatrix}
\begin{pmatrix}H_x\\H_y\end{pmatrix}.
}
$$

矩阵乘法按“每一行乘以输入列，然后相加”计算。第一行结果就是 $Z_{xx}H_x+Z_{xy}H_y$。下标第一位标记**输出电场方向**，第二位标记**输入磁场方向**。这些下标从来都不是测站位置。

### 三、给四个系数随便设数，算一遍

先拿实数练习：

$$
\mathbf Z=\begin{pmatrix}1&2\\3&4\end{pmatrix},
\quad
\mathbf H=\begin{pmatrix}5\\6\end{pmatrix}.
$$

第一行 $E_x=1\times5+2\times6=17$；第二行 $E_y=3\times5+4\times6=39$，所以

$$
\mathbf E=\begin{pmatrix}17\\39\end{pmatrix}.
$$

真实阻抗的四项一般是复数，但计算规则完全一样。**一个 $2\times2$ 复矩阵含有四个复数，也就是八个实数**。之后的维度分析，就是研究八个数之间是否存在特殊关系。

### 四、磁场用 $H$ 还是 $B$？

如果数据文件里记录的是磁感应强度 $\mathbf B$，且采用 $\mathbf B=\mu_0\mathbf H$，那么将 $\mathbf E=\mathbf Z_H\mathbf H$ 改写为

$$
\mathbf E=\underbrace{\frac{\mathbf Z_H}{\mu_0}}_{\mathbf M}\mathbf B.
$$

因此用 $E/H$ 和 $E/B$ 定义的传递系数数值及单位不同。处理软件和教材若采用不同定义，必须先换算后比较。

### 五、第五个通道 $H_z$ 如何使用？

除了水平的四个阻抗元素，还可以估计**垂直磁场传递函数**：

$$
H_z=T_xH_x+T_yH_y.
$$

$T_x,T_y$ 组成 tipper（倾子、垂直传递函数）向量。它反映电磁感应响应的方向性和某些横向异常；它不是把地下电阻率“直接量出来”的工具。绘制 induction arrows（感应箭头）时，还要确认所用箭头约定，有些约定方向相反。

### 六、先认识特殊的一维模型

对水平层状、各向同性且场源可以用平面波近似的模型：

$$
\mathbf Z=
\begin{pmatrix}0&Z\\-Z&0\end{pmatrix}.
$$

这时 $Z_{xx}=Z_{yy}=0$，且 $Z_{xy}=-Z_{yx}$。水平各方向等价，响应具有特定的对称性，对角元因而消失。

**自检：** 如果 $\mathbf Z=\begin{pmatrix}1&2\\3&4\end{pmatrix}$、$\mathbf H=(0,5)^T$，输出是多少？答：$(10,20)^T$。记住只有第二列对结果起作用。

[下一篇：一维、二维、三维与 TE/TM](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06b-dimensionality) · [专题目录](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06-analysis)


---

<a id="sec-06b-dimensionality"></a>
## 第3单元｜数据分析 02｜一维、二维、三维怎么判断？TE 和 TM 又是什么？

电磁场始终存在于三维空间；“一维、二维、三维地下”说的是**电阻率在几个空间方向上发生变化**。先把这件事弄清楚，再看阻抗张量会简单得多。

![不同维数的地下模型](/images/posts/mt-course/23-dimension-responses.svg)

### 一、一维模型为什么最容易？

一维地球的电阻率写作 $\rho=\rho(z)$。同一深度的岩石在水平各处完全相同，怎么转动水平坐标系，材料分布都不会改变。在理想各向同性 1D 情形，阻抗矩阵满足：

$$
\mathbf Z_{1D}=\begin{pmatrix}0&Z\\-Z&0\end{pmatrix},\qquad
T_x=T_y=0.
$$

这里 $Z$ 仍可随频率改变。阻抗对角元为零以及两个非对角元大小相等、符号相反，是重要的数学特征。**但观测有噪声时，实际数据不可能逐项严格等于这些理想数值。**

### 二、二维模型究竟多了什么？

假设沿 $x$ 方向无限延伸的盆地或断层，沿 $x$ 不改变，但沿 $y$ 和深度 $z$ 改变：

$$\rho=\rho(y,z).$$

$x$ 便是地质**走向（strike）**。把坐标轴恰好选在走向及其垂直方向后，理想二维模型有

$$
\mathbf Z_{2D}=\begin{pmatrix}0&Z_{xy}\\Z_{yx}&0\end{pmatrix}.
$$

这次 $Z_{xy}$ 与 $-Z_{yx}$ 通常不再相等。原因是沿走向的电流和跨走向的电流，面对导电结构时经历的物理环境不同。

### 三、TE / TM 为什么是两种模式？

仍让地质走向沿 $x$：

**TE（横向电）模式**的水平电场平行走向，主要研究 $(E_x,H_y)$，因而对应 $Z_{xy}$。

**TM（横向磁）模式**的水平磁场平行走向，主要研究 $(E_y,H_x)$，对应 $Z_{yx}$。

从物理上说，TE 中电流倾向沿着构造走向流动；TM 中电流跨越横向电性界面，更容易在某些导电异常边缘形成电荷积累，产生电流性畸变。两种模式对断层与盆地边界的响应因此可能不同。

要留意：有些资料会采用不同的走向轴定义，导致 TE/TM 与 $xy$、$yx$ 的对应交换。理解极化方向，比单记下标更可靠。

### 四、如果电极线没对准地质走向呢？

现实中测站电极方向由布设条件决定，未必恰好沿地质走向。虽然真实结构仍为二维，**未经旋转的测量坐标中，阻抗对角元却可能不是零**。

这是因为我们用一组“歪着的尺子”描述原本具有对称性的物理响应。稍后我们会用旋转矩阵把张量变回走向坐标。不能看见 $Z_{xx}\ne0$ 就宣称“这里一定是三维”。

### 五、三维模型又是什么？

如果导电体在三个空间方向上都有有限范围，或电阻率满足 $\rho=\rho(x,y,z)$，就需要一般的三维模型。在一般坐标下四个阻抗元素都可能非零，tipper 的两个分量也可能具有明显响应。

但一个某频率的测量张量没有唯一的“地下几维”标签。实际观测可能包含三维浅部畸变、区域二维结构，以及更深的一维背景；随着周期改变，数据敏感的深度范围也会变化。

### 六、一个从矩阵观察维度的小练习

假设某频率的无噪声阻抗为

$$
\mathbf Z=
\begin{pmatrix}0&3+2i\\-3-2i&0\end{pmatrix}.
$$

它符合一维各向同性的特征吗？符合，两个对角元为零、非对角元互为负数。

如果改成

$$
\mathbf Z=
\begin{pmatrix}0&3+2i\\-5-i&0\end{pmatrix},
$$

它更接近在走向坐标中的二维特征。不过仅凭一个频率上的形式，**仍不能唯一证明**地下模型的真实维度。

下一篇我们解决一个令人头疼的问题：为什么变换坐标时，张量要在左右两边各乘一次旋转矩阵？

[继续：旋转矩阵从哪里来？](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06c-rotation) · [专题目录](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06-analysis)


---

<a id="sec-06c-rotation"></a>
## 第4单元｜数据分析 03｜为什么旋转阻抗张量，要左右各乘一个矩阵？

终于轮到一个看起来很吓人的公式：

$$
\boxed{\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T}.
$$

这条式子的关键在于：电、磁场的物理箭头都没有动，变化的只是用来描述它们的坐标轴。那么，为什么阻抗矩阵左右要各乘一次旋转矩阵？

### 一、向量本身不变，分量可以改变

![同一个电场在不同坐标系中的分量](/images/posts/mt-course/22-coordinate-rotation.svg)

设一根箭头指向东北方，长度为 $\sqrt2$。若 $x$ 指北、$y$ 指东，它有分量

$$\mathbf v=\begin{pmatrix}1\\1\end{pmatrix}.$$

现在把新 $x'$ 轴转到箭头方向，新的 $y'$ 轴与其垂直；箭头当然还是原来的箭头，但它在新坐标里的分量会成为 $(\sqrt2,0)$。

所以**坐标变换改变的是分量数字，不是物理场本身**。

### 二、旋转矩阵从两个投影得到

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

### 三、什么叫转置？为什么旋转的逆矩阵特别简单？

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

### 四、真正推导阻抗张量旋转公式

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

### 五、算一个 90° 的具体例子

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

### 六、最后认识三个旋转不变量

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

### 七、自检

1. $\mathbf v=(1,0)^T$ 转动坐标 $90^\circ$ 后是什么？答：$(0,-1)^T$，在本节约定下。
2. 为什么 $\mathbf Z'$ 不是简单的 $\mathbf R\mathbf Z$？因为磁场输入也要变换坐标。
3. 迹与反迹相同是否足以证明两个阻抗张量物理相同？不够，它们只包含部分信息。

[下一节：WAL 不变量与维度判断](/posts/2026/10/mt-course-06d-invariants/) · [专题目录](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06-analysis)
