---
title: "MT 讲义 02｜从麦克斯韦方程一步步推到阻抗"
date: 2026-10-21
permalink: /posts/2026/10/mt-course-02-maxwell-forward/
excerpt: "把旋度、偏导数、简谐场、控制方程、半空间解、视电阻率和层状阻抗递推连成一条可追踪的推导。"
tags: [大地电磁, 麦克斯韦方程, 正演, 阻抗]
---
{% include mt-course-math.html %}

## 1. 先弄清微分符号

普通导数 $dE/dz$ 描述**只有深度一个自变量**时，电场随深度变化率。偏导数 $\partial E/\partial z$ 则表示在电场也依赖 $x,y,t$ 时，先把其余变量固定，再改变 $z$。例如

$$
E(x,z,t)=xz^2\cos(\omega t)
\ \Longrightarrow\
\frac{\partial E}{\partial z}=2xz\cos(\omega t).
$$

$\nabla$ 读作 nabla；$\nabla\cdot\mathbf E$ 是散度，$\nabla\times\mathbf E$ 是旋度；$\nabla^2 E$ 是拉普拉斯算子。对一个标量分量有

$$
\nabla^2E=\frac{\partial^2 E}{\partial x^2}+
\frac{\partial^2 E}{\partial y^2}+
\frac{\partial^2 E}{\partial z^2}.
$$

**理解旋度**：它告诉我们矢量场在一个很小的闭合回路附近，有多强的环流倾向。譬如其 $y$ 分量包含 $\partial E_x/\partial z-\partial E_z/\partial x$。不要把旋度误解成“电场箭头本身每时每刻在转圈”。

## 2. 麦克斯韦方程、材料关系与近似

时域四个麦克斯韦方程是

$$
\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t},\qquad
\nabla\times\mathbf H=\mathbf J+\frac{\partial\mathbf D}{\partial t},
$$

$$
\nabla\cdot\mathbf D=\rho_f,\qquad
\nabla\cdot\mathbf B=0.
$$

此处 $\rho_f$ 表示自由电荷密度，**不要与后面的电阻率 $\rho$ 混淆**。物质本构关系为

$$
\mathbf D=\varepsilon\mathbf E,\quad
\mathbf B=\mu\mathbf H,\quad
\mathbf J=\sigma\mathbf E,\quad \rho=1/\sigma.
$$

$\varepsilon$ 为介电常数，$\mu$ 为磁导率，$\sigma$ 为电导率。对 MT 低频、普通导电大地，经常在介质中假设 $\sigma\gg\omega\varepsilon$，忽略位移电流（具体介质与频率需要检查）。于是

$$
\nabla\times\mathbf H\simeq\sigma\mathbf E.
$$

采用复数简谐场的 $e^{i\omega t}$ 约定，$\partial/\partial t$ 作用在相量上等效于 $i\omega$，令 $\mu\simeq\mu_0$，得到

$$
\nabla\times\mathbf E=-i\omega\mu_0\mathbf H,\qquad
\nabla\times\mathbf H=\sigma\mathbf E.
$$

## 3. 从两个一阶方程得到控制方程

对第一个方程两边取旋度：

$$
\nabla\times(\nabla\times\mathbf E)
=-i\omega\mu_0(\nabla\times\mathbf H)
=-i\omega\mu_0\sigma\mathbf E.
$$

使用矢量恒等式

$$
\nabla\times(\nabla\times\mathbf E)
=\nabla(\nabla\cdot\mathbf E)-\nabla^2\mathbf E.
$$

若讨论**均匀、各向同性且无自由电荷**的区域，$\nabla\cdot\mathbf E=0$，得

$$
-\nabla^2\mathbf E=-i\omega\mu_0\sigma\mathbf E
\quad\Longrightarrow\quad
\boxed{\nabla^2\mathbf E-i\omega\mu_0\sigma\mathbf E=0}.
$$

这是本课程的重要电场控制方程。在电导率变化的区域，不能毫无条件地直接套用上述均匀介质形式；应从原始旋度方程及界面条件处理。对常数 $\sigma,\mu$ 可类似推导磁场方程。

## 4. 为什么“一维”会突然容易许多？

一维层状模型指地层在水平向无限延伸，电阻率仅随 $z$ 改变，因此 $\partial/\partial x=\partial/\partial y=0$。若取一组极化 $(E_x,H_y)$，则

$$
\frac{dE_x}{dz}=-i\omega\mu_0 H_y,\qquad
-\frac{dH_y}{dz}=\sigma E_x.
$$

对第一式再对 $z$ 求导：

$$
\frac{d^2E_x}{dz^2}
=-i\omega\mu_0\frac{dH_y}{dz}
=i\omega\mu_0\sigma E_x.
$$

因此

$$
\boxed{\frac{d^2E_x}{dz^2}-k^2E_x=0,\quad
k^2=i\omega\mu_0\sigma}.
$$

注意这里 $k$ 是**复数传播常数**，与无耗波动中的纯实数波数有联系但并不相同。

## 5. 解这个二阶微分方程

因为指数函数求导后仍是自己乘上常数，所以尝试 $E_x=e^{rz}$；代入方程得 $r^2=k^2$，因此 $r=\pm k$。线性方程的通解为

$$
E_x(z)=Ae^{-kz}+Be^{kz}.
$$

对于向下无限延伸的均匀半空间，选取 $\operatorname{Re}k>0$，排除随 $z\to\infty$ 指数增长的分支，于是 $B=0$：

$$
E_x(z)=Ae^{-kz},\qquad
H_y(z)=-\frac1{i\omega\mu_0}\frac{dE_x}{dz}
=\frac{k}{i\omega\mu_0}E_x.
$$

地表阻抗便是

$$
\boxed{Z_{xy}=\frac{E_x}{H_y}=\frac{i\omega\mu_0}{k}
=\sqrt{i\omega\mu_0\rho}}.
$$

这个公式现在有了完整的来历，而不是凭空背下来的等式。

## 6. 为什么是 45°，又为什么会衰减？

利用 $\sqrt{i}=(1+i)/\sqrt2$，可写

$$
Z_{xy}=\sqrt{\omega\mu_0\rho}\,e^{i\pi/4},\qquad
\arg Z_{xy}=45^\circ.
$$

这是**理想均匀半空间及本讲符号约定**的结论。传播常数可分解为

$$
k=\sqrt{i\omega\mu_0\sigma}=
(1+i)\sqrt{\frac{\omega\mu_0\sigma}{2}}
=\alpha+i\beta,
$$

其中 $\alpha=\beta=\sqrt{\omega\mu_0\sigma/2}$。则

$$
E_x(z,t)=\operatorname{Re}\left[Ae^{i\omega t}e^{-(\alpha+i\beta)z}\right]
=|A|e^{-\alpha z}\cos(\omega t-\beta z+\phi_A).
$$

看清两个因子：$e^{-\alpha z}$ 控制振幅衰减，$-\beta z$ 控制相位随深度变化。趋肤深度 $\delta=1/\alpha$，因此

$$
\delta=\sqrt{\frac{2}{\omega\mu_0\sigma}}.
$$

同一套控制方程就自然解释了频率—探测深度关系。

## 7. 视电阻率与相位

由半空间阻抗模长可得 $|Z|^2=\omega\mu_0\rho$。真实地下往往并非均匀半空间，仍习惯把任意频率的阻抗转成“等效半空间”的量：

$$
\boxed{\rho_a(\omega)=\frac{|Z(\omega)|^2}{\omega\mu_0}},
\qquad\phi(\omega)=\arg Z(\omega).
$$

$\rho_a$ 是**视电阻率**；通常不能把它看作某一确定深度的岩石真电阻率。短周期往往对浅部更敏感，长周期往往对深部更敏感，真实敏感范围取决于介质与采样。

## 8. 层状介质阻抗：从下往上递推

![三层地球](/images/posts/mt-course/04-1d-layers.svg)

一层内既有向下衰减项，又可能有来自界面的等效反射项：

$$
E_j(z)=A_j e^{-k_j(z-z_j)}+B_j e^{k_j(z-z_j)}.
$$

在每层界面要求切向电、磁场连续（无外加面电流），于是可以把下一层的输入阻抗传到上一层。令该层厚度为 $h_j$，本征阻抗 $Z_{0j}=\sqrt{i\omega\mu_0\rho_j}$，$k_j=\sqrt{i\omega\mu_0/\rho_j}$，最底半空间的阻抗就是其本征阻抗；上层的递推公式为

$$
\boxed{
Z_j=Z_{0j}\,
\frac{Z_{j+1}+Z_{0j}\tanh(k_jh_j)}
{Z_{0j}+Z_{j+1}\tanh(k_jh_j)}
}.
$$

这里 $\tanh$ 是双曲正切函数，定义为 $\tanh x=(e^x-e^{-x})/(e^x+e^{-x})$，复数时同样成立。这是层状地球解析正演的核心方法：已知每层 $\rho_j,h_j$，算地表 $Z(\omega)$，再转为 $\rho_a$ 与 $\phi$。它也是检验后面数值正演结果的一个可靠基准。

### 小练习

1. 试把 $\frac{dE_x}{dz}=-i\omega\mu_0H_y$ 改写成 $H_y$ 的表达式。答：$H_y=-\frac{1}{i\omega\mu_0}\frac{dE_x}{dz}$。
2. 为什么必须排除 $Be^{kz}$？答：选取物理衰减分支后，它会在无限深处发散。
3. 为什么 $\rho_a$ 不是特定深度的真实 $\rho$？因为观测阻抗是整个地下电性分布综合作用的结果。

**课件定位：** Ch.2 Part I 第 4–35 页，重点 7–20、24–29 页。**下一讲：**[有限差分与多维正演](/posts/2026/10/mt-course-03-finite-difference/)。