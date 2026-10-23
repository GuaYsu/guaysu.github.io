---
title: "10｜不变量、GB 畸变与相位张量｜合并完整版"
date: 2026-10-21
modified: 2026-10-23
permalink: /posts/2026/10/mt-course-06b-tensor-interpretation/
excerpt: "合并关联内容，公式推导、例题、图片与自检题均完整保留。"
tags: [大地电磁, 数学基础, 自学教程]
---
{% include mt-course-math.html %}

## 本篇目录
- [数据分析 04｜WAL 不变量：不旋转坐标也能看出什么？](#sec-06d-invariants)
- [数据分析 05｜近地表导体怎样扭曲电场：Groom–Bailey 畸变](#sec-06e-distortion)
- [数据分析 06｜相位张量怎么来的：先用数字理解矩阵求逆](#sec-06f-phase-tensor)


---

<a id="sec-06d-invariants"></a>
## 第1单元｜数据分析 04｜WAL 不变量：不旋转坐标也能看出什么？

一份阻抗张量包含四个复数，也就是八个实数。把测量坐标转动之后，这些数通常会变。但如果我们希望判断地下结构是不是二维，最好不要让判断结果只取决于电极线恰好指向哪里。

因此需要**旋转不变量**。我们先从最简单的迹、反迹开始，再进入 WAL（Weaver–Agarwal–Lilley）不变量。

### 一、什么是与旋转方向无关的量？

上一篇已经推导

$$\mathbf Z'=\mathbf R\mathbf Z\mathbf R^T.$$

虽然 $Z_{xx}$ 等单个元素会变，但

$$
T=Z_{xx}+Z_{yy},\qquad A=Z_{xy}-Z_{yx}
$$

在纯水平旋转下保持不变。这是一个很实用的思路：将好几个元素组合，消去坐标角度带来的变化。

别忘了复数需要同时考虑实部和虚部。例如 $Z_{xy}=3+2i$，其中实部为 3、虚部为 2。四个复数合起来是八个实数。

### 二、把八个实数重新排成四个复数组合

定义

$$
\begin{aligned}
c_1&=\frac12(Z_{xx}+Z_{yy})=\xi_1+i\eta_1,\\
c_2&=\frac12(Z_{xy}+Z_{yx})=\xi_2+i\eta_2,\\
c_3&=\frac12(Z_{xx}-Z_{yy})=\xi_3+i\eta_3,\\
c_4&=\frac12(Z_{xy}-Z_{yx})=\xi_4+i\eta_4.
\end{aligned}
$$

这里 $\xi$ 读作 xi，$\eta$ 读作 eta，分别是对应组合的**实部与虚部**。八个实数仍然来自原先的四个复阻抗元素，只是换了一种排列方式。

从定义一眼看出 $c_1$ 是迹的一半，$c_4$ 是反迹的一半。其余组合包含对角元之差和非对角元之和。

### 三、WAL 前四个不变量怎么计算？

令

$$
\boxed{I_1=\sqrt{\xi_1^2+\xi_4^2}},\qquad
\boxed{I_2=\sqrt{\eta_1^2+\eta_4^2}}.
$$

它们提供实部、虚部的重要尺度，但本身并不能独立判定地下维度。

定义

$$
\boxed{I_3=\frac{\sqrt{\xi_2^2+\xi_3^2}}{I_1}},\qquad
\boxed{I_4=\frac{\sqrt{\eta_2^2+\eta_3^2}}{I_2}}.
$$

$I_3,I_4$ 表示偏离理想一维形式的某种程度：对于无噪声的一维各向同性响应，$c_2=c_3=0$，因此 $I_3=I_4=0$，前提是对应分母非零。

为什么会出现这些组合？在旋转坐标时，$(\xi_2,\xi_3)$ 会像一个二维向量一样混合，但它们长度 $\sqrt{\xi_2^2+\xi_3^2}$ 不变；$(\eta_2,\eta_3)$ 同理。因此构造这种模长比值，就能减少对电极方向的依赖。

### 四、后面三个不变量为何复杂？

为进一步检测复数响应中实部和虚部的相对方向、以及电流性畸变相关特征，先构造两组归一化组合。

对 $i,j\in\{1,2,3,4\}$，令

$$
d_{ij}=
\frac{\xi_i\eta_j-\xi_j\eta_i}{I_1I_2},
\qquad
s_{ij}=
\frac{\xi_i\eta_j+\xi_j\eta_i}{I_1I_2}.
$$

在某些写法中，分母也会显式写成 $\sqrt{\xi_1^2+\xi_4^2}\sqrt{\eta_1^2+\eta_4^2}$，意思相同。数学上分子像二维向量行列式与对称组合，可以反映两个方向之间的相对关系。

在本节使用的 WAL 约定下：

$$
I_5=s_{41},\qquad I_6=d_{41},
$$

$$
\boxed{
I_7=\frac{d_{41}-d_{23}}
{\sqrt{(d_{12}-d_{34})^2+(d_{13}+d_{24})^2}}.
}
$$

此外可以定义派生指标

$$
Q=\sqrt{(d_{12}-d_{34})^2+(d_{13}+d_{24})^2}.
$$

这样复杂的一串式子为什么有用？它将原本八个依赖坐标的实数提取成七个旋转不变量，剩余一个自由度可以对应方向角。**计算时一定要注意分母为零或数值很小时，某些指标会无定义或非常不稳定**。不能把无定义的结果机械记成 0。

### 五、这些指标在一维、二维时有什么特征？

对于理想、无噪声的一维介质：$I_3$ 到 $I_6$ 为零，$I_7$ 可能无定义；对理想二维介质，$I_5$、$I_6$ 为零，$I_7$ 为零或无定义。非理想畸变和实际噪声会破坏这些严格等式。

不同研究也曾使用 $0.1$ 或 $0.15$ 之类的经验阈值，将一个小数值视为“统计上接近零”。这些阈值与数据不确定度、处理方法相关，不应视为普适的物理常数。

### 六、经典 Swift skew 又是什么？

一个计算很容易的经典指标是

$$
\boxed{
\mathrm{skew}_{S}=
\frac{|Z_{xx}+Z_{yy}|}{|Z_{xy}-Z_{yx}|}.
}
$$

例如 $Z_{xx}=0.1$、$Z_{yy}=0.1$、$Z_{xy}=4$、$Z_{yx}=-4$，则

$$
\mathrm{skew}_{S}=\frac{0.2}{8}=0.025.
$$

它很小，与某种近一维或二维的简洁响应相容；不过浅部畸变和噪声也会影响它。经典经验上曾使用 0.2 作为一个粗略门槛，但不能因为某点低于 0.2 就宣布该地区一定二维。

更复杂的 Bahr phase-sensitive skew 试图利用相位信息减少某类畸变影响；仍需要注意噪声较高时的稳定性。WAL 系列提供更系统的组合指标，但它们也无法单独解决全部地质结构判断。

### 七、读完需要真正记住什么？

不需要一口气背下 $I_1$ 到 $I_7$ 的所有下标。请先理解这三层逻辑：

1. 四个复阻抗元素是原始输入；
2. 将它们分成迹、反迹及对称/反对称组合，构造随旋转保持不变的量；
3. 再结合这些量及其误差，判断哪些简化模型有希望解释数据。

实际做维度分析时，还要同时考虑相位张量、区域地质走向、tipper、频率依赖与不同测点一致性。

[下一节：Groom–Bailey 畸变分解](/posts/2026/10/mt-course-06b-tensor-interpretation/#sec-06e-distortion) · [专题目录](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06-analysis)


---

<a id="sec-06e-distortion"></a>
## 第2单元｜数据分析 05｜近地表导体怎样扭曲电场：Groom–Bailey 畸变

想象我们真正关心的是十几千米深处的大型导电结构。可测站正好落在一小块浅部异常导体旁边，这个小体的尺度远小于区域构造，却会改变局部电流的方向和强弱。

于是仪器测到的电场，可能是**区域电场经过局部结构重新分配以后形成的结果**。

![浅部导体与局部电流畸变](/images/posts/mt-course/24-galvanic-distortion.svg)

### 一、为什么“近地表异常”可能产生电场畸变？

电导率不同的两个材料相接时，电流连续性和电磁边界条件可能导致界面上积累电荷，从而产生额外的局部电场。若异常体尺度远小于研究频带的特征电磁变化尺度，局部响应中的一部分可以近似看成频率无关的**电流性畸变**（galvanic distortion）。

在这种理想近似里：

$$
\boxed{\mathbf E_{\rm obs}=\mathbf C\mathbf E_{\rm reg}}.
$$

$\mathbf E_{\rm reg}$ 是我们希望知道的区域电场；$\mathbf E_{\rm obs}$ 是测到的电场；$\mathbf C$ 是 $2\times2$ **实数矩阵**。

举个简单例子：

$$
\mathbf C=\begin{pmatrix}2&0\\0&1\end{pmatrix}.
$$

若区域电场是 $(3,4)^T$，那么观测电场就成为 $(6,4)^T$。你看，原本相同的物理区域场，在测点上被不同方向的比例放大了。

### 二、为什么阻抗也会被畸变？

如果暂时假设磁场不发生明显畸变，且区域电场满足

$$
\mathbf E_{\rm reg}=\mathbf Z_{\rm reg}\mathbf H,
$$

代入上一式：

$$
\mathbf E_{\rm obs}
=\mathbf C\mathbf Z_{\rm reg}\mathbf H.
$$

因此测得的阻抗是

$$
\boxed{\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg}}.
$$

因此即使区域结构接近二维，实测阻抗仍可能出现明显的对角项。浅部导电小体会使水平电场的两个方向发生混合，相当于让两个通道互相‘串门’。

### 三、GB 分解试图怎样认识矩阵 C？

Groom–Bailey（GB）分解把局部畸变视为若干具有较直观作用的运算组合。常见符号写作

$$
\boxed{\mathbf C=g\,\mathbf T\,\mathbf S\,\mathbf A},
$$

其中 $g$ 为 site gain（测站增益），$\mathbf T$ 为 twist（扭转），$\mathbf S$ 为 shear（剪切），$\mathbf A$ 为 anisotropy（各向异性增益）。不同文献使用的规范化和矩阵排列约定需要核对。

为了先理解数学，取一种常用的示意参数化。扭转可写作一个旋转类矩阵：

$$
\mathbf T(t)=\frac1{\sqrt{1+t^2}}
\begin{pmatrix}1&-t\\t&1\end{pmatrix}.
$$

它表示两个水平电场分量因局部效应而发生类似旋转的混合。参数 $t$ 控制这种“扭”的程度；其对应角度可由反正切表示。

剪切可以示意为

$$
\mathbf S(s)=\frac1{\sqrt{1+s^2}}
\begin{pmatrix}1&s\\s&1\end{pmatrix}.
$$

于是如果原区域电场为 $(E_x,E_y)^T$，剪切后的两个分量为

$$
E_x'=\frac{E_x+sE_y}{\sqrt{1+s^2}},\qquad
E_y'=\frac{sE_x+E_y}{\sqrt{1+s^2}}.
$$

两条式子告诉我们：**测到的每个方向，都混入了原本另一个方向的场**。

各向异性增益可以用类似

$$
\mathbf A(a)=\frac1{\sqrt{1-a^2}}
\begin{pmatrix}1+a&0\\0&1-a\end{pmatrix},\qquad |a|<1
$$

的形式表示，它使两个方向的幅值按不同的比例变化。标量 $g$ 则使整个电场统一放大或缩小。上述因子可以采用不同的归一化形式，部分比例会被吸收入 $g$。

### 四、为什么把四个效应分开仍很难唯一求出？

假设只在一个频率、一个测站测到了四个复阻抗元素，共八个实数。我们不仅想恢复区域二维阻抗的两个复数分量，还要估计区域走向、局部扭转、剪切以及增益等参数。

未知量过多，而且某些因子可互相补偿。例如整体放大 $g$ 与区域阻抗整体模长的改变，在单站数据里就可能无法区分。这意味着**测到一个完整阻抗张量，并不保证我们能唯一恢复全部畸变参数**。

常见策略是跨多个频率与多个测点联合拟合，并对局部畸变的频率无关性、区域二维走向的一致性等提出约束。即便拟合效果好，绝对 site gain 通常仍可能存在不可辨识问题。

### 五、畸变不是一种单一现象

电流性（galvanic）畸变侧重界面电荷及局部电流重新分配；电感性（inductive）畸变与频率相关的感应电流、次生磁场有关。

同一个浅部导体对不同频率和 TE/TM 模式可能产生不同影响。将所有异常响应简单叫作“静态位移”，会混淆幅值变化、方向混合与频率相关的效应。

**自检：** 若 $\mathbf C=\begin{pmatrix}1&0.5\\0&1\end{pmatrix}$，区域场为 $(2,4)^T$，观测场是多少？答：$(4,4)^T$。因为 $E_x$ 混进了 $0.5E_y$。

[下一节：相位张量为何能抵消实数电流性畸变？](/posts/2026/10/mt-course-06b-tensor-interpretation/#sec-06f-phase-tensor) · [专题目录](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06-analysis)


---

<a id="sec-06f-phase-tensor"></a>
## 第3单元｜数据分析 06｜相位张量怎么来的：先用数字理解矩阵求逆

我们在上篇见到一个麻烦：观测到的阻抗可能是

$$
\mathbf Z_{\rm obs}=\mathbf C\mathbf Z_{\rm reg},
$$

其中 $\mathbf C$ 表示近地表的实数电流性畸变矩阵。能不能构造某个数学量，利用矩阵运算让 $\mathbf C$ 自动抵消？

这正是**相位张量（phase tensor）**的重要动机。为了真正理解它，我们先从最容易算的二阶矩阵开始。

### 一、为什么要把复数矩阵分开？

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

### 二、矩阵的“倒数”究竟是什么？

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

### 三、亲手算一个相位张量

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

### 四、最重要的证明：为什么畸变会抵消？

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

逆矩阵的顺序可以直接验算：$(\mathbf A\mathbf B)(\mathbf B^{-1}\mathbf A^{-1})=\mathbf A\mathbf I\mathbf A^{-1}=\mathbf I$。因此乘积的逆要写成 $\mathbf B^{-1}\mathbf A^{-1}$。矩阵乘法有先后顺序，这一点和普通数字乘法不同。

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

### 五、相位张量元素怎么实际计算？

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

### 六、为什么要把相位张量画成椭圆？

![相位张量的圆和椭圆示意](/images/posts/mt-course/25-phase-tensor-ellipse.svg)

一个二维矩阵可以描述不同方向上的变换。对相位张量常采用主轴、椭圆率和偏斜角绘制椭圆：椭圆长短轴反映不同主方向的响应差异，方向角描述其方位，颜色还可表示偏斜或其他指标。

理想一维各向同性结构的相位张量与单位矩阵成比例，通常绘成圆。理想二维模型在走向坐标中无偏斜，主轴方向与地质走向及其垂直方向有关。对于实际资料，椭圆拉长还可能反映噪声、各向异性或三维结构，因此要结合偏斜角、tipper 和邻近测站一起判断。

### 七、四个值得认识的数

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

### 八、相位张量不能做什么？

相位张量只有四个实数分量，原始复阻抗则包含八个实数分量。构造相位张量时，部分幅度信息已被消去；它适合诊断电性维数与方向性，恢复完整阻抗或进行幅相联合反演时仍需要原始数据。

**练习：** 若 $\mathbf X=2\mathbf I$、$\mathbf Y=3\mathbf I$，则 $\boldsymbol\Phi=(3/2)\mathbf I$，对应理想的各向同性相位张量形态。

[下一节：走向估计与维度诊断](/posts/2026/10/mt-course-06c-strike-static/#sec-06g-strike) · [专题目录](/posts/2026/10/mt-course-06a-tensor-dimensions/#sec-06-analysis)
