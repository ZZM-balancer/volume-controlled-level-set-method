# 弹性结构第一特征值水平集优化

本文档把数学模型说明和 FreeFEM 代码融合在一起。所有 `~~~freefem` 代码块从上到下拼接后，就是一个完整的 `.edp` 程序。

本程序在你原来的水平集拓扑优化框架上，把状态方程由静力弹性方程改为弹性结构广义特征值问题。目标由 compliance 最小化改为第一特征值最大化。速度场仍然通过保体积约束的 Riesz 表示求得，因此不需要在目标函数中另外人为加入任意体积乘子。

## 1. 必要的 load

这里只保留特征值优化和水平集更新需要的模块。

- `distance` 用于重新初始化水平集函数；
- `iovtk` 用于输出 `.vtu` 文件，方便 ParaView 后处理；
- 不使用 `PETSc`，矩阵求解全部使用 FreeFEM 原生 `matrix`。

~~~freefem
load "distance"
load "iovtk"
verbosity = 0;
~~~

## 2. 参数设置

| 变量名 | 含义 |
|:---:|:---|
| `Titer` | 每次形状更新的伪时间步长 |
| `iterMax` | 最大优化迭代步数 |
| `mu` | 线弹性 Lamé 常数 \(\mu\) |
| `lameLambda` | 线弹性 Lamé 常数 \(\lambda_L\)。这里不用 `lambda` 作为变量名，避免与特征值概念混淆 |
| `weak` | 刚度矩阵中的弱材料系数，即 \(\rho_K=\mathrm{weak}+(1-\mathrm{weak})H_\gamma(\phi)\) |
| `weakM` | 质量矩阵中的弱材料系数。这里取 `1e-4` |
| `blockMass` | 右端小黑块在质量矩阵中的等效密度，用于模拟端部附加质量 |
| `gamma` | 光滑 Heaviside 函数的界面宽度 |
| `sigma` | 特征值求解的 shift 参数。这里取 \(0.0\)，用于寻找靠近 \(0.0\) 的最小正特征值 |
| `nev` | 计算的特征值个数。虽然当前只优化第一个特征值，但多算几个有利于观察模态交换 |
| `object` | 优化过程信息输出文件 |
| `ModelMatlabName` | 输出给 Matlab 画图的数据文件 |

目标是最大化第一特征值 \(\lambda_1\)。实际速度场方程按最小化

\[
J(\phi)=-\lambda_1(\phi)
\]

来写，因此速度场方向是使 \(J\) 下降，也就是使 \(\lambda_1\) 上升的方向。

~~~freefem
real Titer = 0.06;
int iterMax = 900;

real mu = 8.0;
real lameLambda = 1.0;
real weak = 1e-2;
real weakM = 1e-4;
real blockMass = 35.0;
real gamma = 0.005;

real sigma = 0.0;
int nev = 4;
real modeGapTol = 1e-2;

ofstream object("./Output/MainEigen.obj");
string ModelMatlabName = "matlab/Model1EigenShape1.m";
~~~

## 3. 水平集函数和密度函数

设计区域为固定区域 \(D\)。材料区域用水平集函数表示：

\[
\Omega=\{x\in D:\phi(x)<0\}.
\]

为了在固定区域上计算弹性特征值问题，这里区分刚度密度和质量密度。刚度矩阵使用

\[
\rho_K(\phi)=\mathrm{weak}+(1-\mathrm{weak})H_\gamma(\phi),
\]

质量矩阵使用

\[
\rho_M(\phi)=\mathrm{weakM}+(1-\mathrm{weakM})H_\gamma(\phi).
\]

本程序按你的要求取

\[
\mathrm{weakM}=1e-4,
\]

也就是空洞区域不提供质量；而刚度矩阵仍保留很小的 ersatz 刚度，用于避免刚度矩阵奇异。\(\phi<0\) 对应真实材料区域，\(\phi>0\) 对应空洞或弱材料区域。

采用的刚度密度光滑 Heaviside 函数为

\[
H_{\gamma,\mathrm{weak}}(\phi)
=
\begin{cases}
\mathrm{weak}, & \phi>\gamma,\\[1mm]
(1-\mathrm{weak})\left(\dfrac12-\dfrac{\phi}{\gamma}+\dfrac{\phi^2}{2\gamma^2}\right)+\mathrm{weak},
& 0<\phi\le \gamma,\\[2mm]
(1-\mathrm{weak})\left(\dfrac12-\dfrac{\phi}{\gamma}-\dfrac{\phi^2}{2\gamma^2}\right)+\mathrm{weak},
& -\gamma\le \phi\le 0,\\[2mm]
1, & \phi< -\gamma.
\end{cases}
\]

~~~freefem
macro HeaviSide(phi)
(
    (phi > gamma) * weak
  + (0 < phi && phi <= gamma)
    * ((1.0 - weak) * (0.5 - phi/gamma + 0.5*phi^2/gamma^2) + weak)
  + (-gamma <= phi && phi <= 0)
    * ((1.0 - weak) * (0.5 - phi/gamma - 0.5*phi^2/gamma^2) + weak)
  + (phi < -gamma) * 1.0
) //

macro HeaviSideMass(phi)
(
    (phi > gamma) * weakM
  + (0 < phi && phi <= gamma)
    * ((1.0 - weakM) * (0.5 - phi/gamma + 0.5*phi^2/gamma^2) + weakM)
  + (-gamma <= phi && phi <= 0)
    * ((1.0 - weakM) * (0.5 - phi/gamma - 0.5*phi^2/gamma^2) + weakM)
  + (phi < -gamma) * 1.0
) //

macro Hprime(phi)
(
    (0 < phi && phi <= gamma)
    * ((1.0 - weak) * (-1.0/gamma + phi/(gamma^2)))
  + (-gamma <= phi && phi <= 0)
    * ((1.0 - weak) * (-1.0/gamma - phi/(gamma^2)))
) //

macro BlackBlock()
(
    (x >= 0.9 && x <= 1.0 && y >= -0.05 && y <= 0.05)
) //

macro RhoWithBlock(phi)
(
    BlackBlock * 1.0
  + (1.0 - BlackBlock) * HeaviSide(phi)
) //

macro RhoMassWithBlock(phi)
(
    BlackBlock * blockMass
  + (1.0 - BlackBlock) * HeaviSideMass(phi)
) //

~~~

这里额外加入了一个贴在右端中部边界 `bdrightNN` 上的小黑色块区域

\[
D_{\mathrm{block}}=[0.9,1]\times[-0.05,0.05].
\]

该区域作为非设计区域处理，但刚度和质量的含义不同。刚度矩阵中，小黑块仍是普通实心材料，

\[
\rho_K(x)=1,\qquad x\in D_{\mathrm{block}}.
\]

质量矩阵中，小黑块用于模拟端部附加质量，因此取

\[
\rho_M(x)=\mathrm{blockMass},\qquad x\in D_{\mathrm{block}},
\]

本程序默认 `blockMass = 10.0`。在其余可设计区域中，刚度密度取 \(H_{\gamma,\mathrm{weak}}(\phi)\)，质量密度取 \(H_{\gamma,\mathrm{weakM}}(\phi)\)，且这里 \(\mathrm{weakM}=0\)。因此实际使用

\[
\rho_K(\phi)=\chi_{D_{\mathrm{block}}}+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}(\phi),
\]

\[
\rho_M(\phi)=\mathrm{blockMass}\,\chi_{D_{\mathrm{block}}}+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weakM}}(\phi).
\]

体积监测和保体积约束仍使用刚度/实体密度 \(\rho_K\)，不会因为 `blockMass` 上调而改变初始体积。速度场和体积约束只对设计区域起作用，所以形状导数中使用

\[
\nabla \rho_{K,\mathrm{design}}
=(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}'(\phi)\nabla\phi.
\]

## 4. 线弹性记号

位移场为

\[
u=(u_x,u_y).
\]

线弹性应变张量为

\[
\varepsilon(u)=\frac12(\nabla u+\nabla u^T).
\]

在二维情形下，代码中使用 Voigt 形式

\[
\varepsilon(u)
=
\left(
\partial_x u_x,
\partial_y u_y,
\frac{\partial_y u_x+\partial_x u_y}{\sqrt2}
\right).
\]

这样有

\[
\varepsilon(u):\varepsilon(v)
=
\left[
\partial_x u_x,
\partial_y u_y,
\frac{\partial_y u_x+\partial_x u_y}{\sqrt2}
\right]
\cdot
\left[
\partial_x v_x,
\partial_y v_y,
\frac{\partial_y v_x+\partial_x v_y}{\sqrt2}
\right].
\]

散度为

\[
\operatorname{div}u=\partial_xu_x+\partial_yu_y.
\]

~~~freefem
macro epsilon(ux, uy) [dx(ux), dy(uy), (dy(ux) + dx(uy))/sqrt(2.0)] //
macro grad(u) [dx(u), dy(u)] //
macro Grad(ux, uy) [dx(ux), dy(ux), dx(uy), dy(uy)] //
macro div(ux, uy) (dx(ux) + dy(uy)) //

macro gradRhoDesign(phi)
(
    (1.0 - BlackBlock) * Hprime(phi) * grad(phi)
) //
~~~

## 5. 网格和有限元空间

计算区域为二维矩形

\[
D=(-1,1)\times(-0.5,0.5).
\]

边界标号沿用你原来的设置：

| label | 用途 |
|:---:|:---|
| `2` | 左边界，弹性位移固定边界 \(\Gamma_D\) |
| `3` | 右侧中间边界，用于速度场边界控制 |
| `11` | 左右方向速度约束边界，用于约束 \(\theta_x\) |
| `12` | 上下方向速度约束边界，用于约束 \(\theta_y\) |

有限元空间选取为：

| 空间 | 有限元 | 作用 |
|:---:|:---:|:---|
| `Uh` | \([P2,P2]\) | 弹性特征振型空间 |
| `Thetah` | \([P1,P1]\) | 水平集速度场空间 |
| `P1h` | \(P1\) | `distance` 重新初始化的输出空间 |
| `Phih` | \(P1dc\) | dual-P^1 DG 水平集函数空间 |

这里要注意：`distance` 可以接受一般有限元函数或者表达式作为输入，但输出数组必须属于 `P1h`。因此重新初始化时先输出到 `P1h phiRe`，再赋值回 `Phih`。

~~~freefem
border bdrightN1(t = -0.5, -0.05) { x=1;  y=t; label=11; };
border bdrightNN(t = -0.05, 0.05) { x=1;  y=t; label=3;  };
border bdrightN2(t = 0.05, 0.5)   { x=1;  y=t; label=11; };
border bdleftD(t = 0.5, -0.5)     { x=-1; y=t; label=2;  };
border bdupN(t = 1, -1)           { x=t;  y=0.5; label=12; };
border bdbottomN(t = -1, 1)       { x=t;  y=-0.5; label=12; };

int n = 4;
mesh Th = buildmesh(
      bdrightN1(9*n) + bdrightNN(3*n) + bdrightN2(9*n)
    + bdleftD(20*n) + bdupN(40*n) + bdbottomN(40*n)
);

savevtk("./figure/Th1.vtu", Th);

fespace Uh(Th, [P2, P2]);
fespace Thetah(Th, [P1, P1]);
fespace P1h(Th, P1);
fespace Phih(Th, P1dc);

P1h h1 = hTriangle;
object << "maximal size of a mesh triangle h = " << h1[].max << endl;
object << "minimal size of a mesh triangle h = " << h1[].min << endl;
~~~

## 6. 初始形状

初始材料区域仍采用多个圆形区域的并集。先构造一个函数 `phi0`，使得

\[
\phi_0<0
\]

对应材料区域。然后通过 `distance` 重新初始化为近似有向距离函数。

需要强调的是，本程序中 `phi` 的主要演化空间是 `Phih=P1dc`，但 `distance` 的输出必须是 `P1h`：

\[
\texttt{distance(Th, phi0, phiInit[])}.
\]

然后再令

\[
\phi=\phi_{\mathrm{init}}.
\]

~~~freefem
macro R1() (0.2 - sqrt((abs(x))^2 + (abs(y) - 0.5)^2)) //
macro R2() (0.2 - sqrt((abs(x) - 0.6)^2 + (abs(y) - 0.5)^2)) //
macro R3() (0.1 - sqrt((abs(x) - 0.3)^2 + (abs(y) - 0.25)^2)) //
macro R4() (0.1 - sqrt((abs(x) - 1.0)^2 + (abs(y) - 0.25)^2)) //
macro R5() (0.2 - sqrt((abs(x))^2 + (abs(y))^2)) //
macro R6() (0.2 - sqrt((abs(x) - 0.6)^2 + (abs(y))^2)) //

func phi0 = (R1 >= 0) * R1
          + (R2 >= 0) * R2
          + (R3 >= 0) * R3
          + (R4 >= 0) * R4
          + (R5 >= 0) * R5
          + (R6 >= 0) * R6
          + (R1 < 0 && R2 < 0 && R3 < 0 && R4 < 0 && R5 < 0 && R6 < 0)
            * max(R1, max(R2, max(R3, max(R4, max(R5, R6)))));

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho, RhoMass;
phi = phiInit;
Rho = RhoWithBlock(phi);
RhoMass = RhoMassWithBlock(phi);

plot(Rho, fill = 1, value = 1, cmm = "Rho0", wait = false);
savevtk("./figure/Rho0.vtu", Th, Rho);
savevtk("./figure/RhoMass0.vtu", Th, RhoMass);
savevtk("./figure/phi0.vtu", Th, phi);
~~~

## 7. 弹性结构第一特征值问题

给定当前水平集函数 \(\phi\)，求广义特征值问题

\[
K(\rho_K)u=\lambda M(\rho_M)u.
\]

弱形式为：求 \((\lambda,u)\)，使得

\[
a_{\rho_K}(u,v)=\lambda m_{\rho_M}(u,v),
\qquad \forall v\in [H^1_{\Gamma_D}(D)]^2.
\]

其中

\[
a_{\rho_K}(u,v)
=
\int_D
\rho_K
\left(
2\mu\varepsilon(u):\varepsilon(v)
+
\lambda_L\operatorname{div}u\operatorname{div}v
\right)dx,
\]

\[
m_{\rho_M}(u,v)
=
\int_D
\rho_M u\cdot v\,dx.
\]

归一化条件为

\[
m_{\rho_M}(u,u)=1.
\]

在代码中，`EigenOP` 实际组装的是 shift 后的矩阵

\[
K(\rho)-\sigma M(\rho).
\]

这里取 `sigma=0`，所以它就是刚度矩阵 \(K(\rho)\)。Dirichlet 边界条件应当在刚度矩阵和质量矩阵中一致处理，即二者都对应同一个满足 \(u=0\) 的离散空间。实际在 FreeFEM 中，为避免边界自由度产生靠近目标特征值的伪特征值，刚度型矩阵和质量矩阵使用不同的 `tgv`：刚度矩阵取很大的 `tgv=1e30`，质量矩阵取 `tgv=1.0`。这样边界自由度对应的伪特征值约为 \(10^{30}\)，不会影响靠近零的第一正特征值。

~~~freefem
Uh [modeX, modeY], [testX, testY];
real lambda1, lambda2, lambda3, lambda4;
real gap12, gap12Ratio;

varf EigenOP([ux, uy], [vx, vy]) =
int2d(Th)(
      Rho * (
          2.0 * mu * (epsilon(ux, uy)' * epsilon(vx, vy))
        + lameLambda * div(ux, uy) * div(vx, vy)
      )
    /*- sigma * RhoMass * (ux * vx + uy * vy)*/
)
+ on(2, ux = 0, uy = 0);

varf Mass([ux, uy], [vx, vy]) =
int2d(Th)(RhoMass * (ux * vx + uy * vy))
+ on(2, ux = 0, uy = 0);

int showeigsize = 1;
macro SolveFirstEigen(modeX, modeY, lambda1)
{
    matrix OPm = EigenOP(Uh, Uh, tgv = 1e30);
    matrix Mm = Mass(Uh, Uh, tgv = 1e10);
    set(OPm, solver = sparsesolver);
    if (showeigsize == 1)
    {
        cout << "EigenOP matrix size = " << OPm.n << " x " << OPm.m << endl;

        object << "EigenOP matrix size = " << OPm.n << " x " << OPm.m << endl;

        showeigsize = 0;
    }

    real[int] eigVal(nev);
    Uh[int] [eigX, eigY](nev);

    int nconv = EigenValue(OPm, Mm,
                           sym = true,
                           sigma = sigma,
                           value = eigVal,
                           vector = eigX,
                           tol = 1e-9,
                           maxit = 0,
                           ncv = 0);

    if (nconv < 1)
    {
        cout << "EigenValue failed: no converged eigenvalue." << endl;
        exit(1);
    }

    lambda1 = eigVal[0];
    lambda2 = 0.0;
    lambda3 = 0.0;
    lambda4 = 0.0;
    gap12 = 1e30;
    gap12Ratio = 1e30;

    real denomGap = abs(lambda1);
    if (denomGap < 1e-30)
    {
        denomGap = 1e-30;
    }
    if (nconv > 1)
    {
        lambda2 = eigVal[1];
        gap12 = lambda2 - lambda1;
        gap12Ratio = gap12 / denomGap;
    }
    if (nconv > 2)
    {
        lambda3 = eigVal[2];
    }
    if (nconv > 3)
    {
        lambda4 = eigVal[3];
    }

    cout << "Eigenvalues: lambda1 = " << lambda1
         << ", lambda2 = " << lambda2
         << ", lambda3 = " << lambda3
         << ", lambda4 = " << lambda4 << endl;
    cout << "gap12 = " << gap12
         << ", gap12/lambda1 = " << gap12Ratio << endl;

    object << "Eigenvalues: lambda1 = " << lambda1
           << ", lambda2 = " << lambda2
           << ", lambda3 = " << lambda3
           << ", lambda4 = " << lambda4 << endl;
    object << "gap12 = " << gap12
           << ", gap12/lambda1 = " << gap12Ratio << endl;

    if ((nconv > 1) && (gap12Ratio < modeGapTol))
    {
        cout << "Warning: possible mode switching or nearly multiple first eigenvalue." << endl;
        object << "Warning: possible mode switching or nearly multiple first eigenvalue." << endl;
    }

    modeX[] = eigX[0][];

    real massNorm = sqrt(int2d(Th)(RhoMass * (modeX^2 + modeY^2)));
    modeX[] = modeX[] / massNorm;

    cout << "EigenOP Success! "<< endl;
} //
~~~

## 8. 第一特征值的形状导数

现在推导第一特征值的导数。假设 \(\lambda_1\) 是简单特征值，对应振型为 \(u_1\)，并满足归一化条件

\[
m_{\rho_M}(u_1,u_1)=1.
\]

特征值方程为

\[
a_{\rho_K}(u_1,v)=\lambda_1 m_{\rho_M}(u_1,v).
\]

对 \(\rho\) 求导，由于特征值问题是自伴随的，振型导数项会通过特征值方程抵消，得到

\[
d\lambda_1
=
a_{\rho_K}'(u_1,u_1)-\lambda_1m_{\rho_M}'(u_1,u_1).
\]

又因为

\[
a_{\rho_K}'(u_1,u_1)
=
\int_D
\left(
2\mu\varepsilon(u_1):\varepsilon(u_1)
+
\lambda_L(\operatorname{div}u_1)^2
\right)
\delta\rho_K\,dx,
\]

\[
m_{\rho_M}'(u_1,u_1)
=
\int_D |u_1|^2\delta\rho_M\,dx,
\]

由于

\[
\delta\rho_K=(1-\mathrm{weak})H_\gamma'(\phi)\delta\phi,
\qquad
\delta\rho_M=(1-\mathrm{weakM})H_\gamma'(\phi)\delta\phi,
\]

把 \(\delta\rho_K\) 提出来，得到

\[
d\lambda_1
=
\int_D
S(u_1,\lambda_1)\delta\rho_K\,dx,
\]

其中

\[
S(u_1,\lambda_1)
=
2\mu\varepsilon(u_1):\varepsilon(u_1)
+
\lambda_L(\operatorname{div}u_1)^2
-
\lambda_1\frac{1-\mathrm{weakM}}{1-\mathrm{weak}}|u_1|^2.
\]

这里按当前代码中的 ersatz 密度归一化方式，把质量项写为 \(\lambda_1(1-\mathrm{weakM})|u_1|^2/(1-\mathrm{weak})\)。

水平集输运方程为

\[
\phi_t+\theta\cdot\nabla\phi=0.
\]

在非设计小黑块外，密度随水平集变化；在小黑块内，密度固定为 \(1\)。因此形状导数中真正参与变化的是

\[
\nabla\rho_{K,\mathrm{design}}
=(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}'(\phi)\nabla\phi.
\]

沿水平集输运方程有

\[
\rho_t
=-(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}'(\phi)\theta\cdot\nabla\phi
=-\theta\cdot\nabla\rho_{K,\mathrm{design}}.
\]

因此

\[
d\lambda_1(\theta)
=-\int_D
S(u_1,\lambda_1)
\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

令

\[
J(\phi)=-\lambda_1(\phi),
\]

则

\[
dJ(\theta)
=
\int_D
S(u_1,\lambda_1)
\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

## 9. 保体积速度场

你这里采用速度场强制体积守恒，所以不需要另外在目标函数中手动加入任意体积乘子。体积为

\[
\operatorname{Vol}(\phi)
=
\int_D \rho(\phi)\,dx.
\]

沿水平集输运方程演化时，小黑块部分固定不变，所以

\[
\frac{d}{dt}\operatorname{Vol}(\phi)
=
\int_D \rho_t\,dx
=-\int_D \theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

因此只要速度场满足

\[
\int_D \theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0,
\]

就得到连续层面的体积守恒。

于是求速度场 \(\theta\) 和约束变量 \(p\)，使得

\[
\int_D \nabla\theta:\nabla w\,dx
+
\int_D \theta\cdot w\,dx
+p\int_D w\cdot\nabla\rho_{K,\mathrm{design}}\,dx
=-dJ(w),
\]

即

\[
\int_D \nabla\theta:\nabla w\,dx
+
\int_D \theta\cdot w\,dx
+p\int_D w\cdot\nabla\rho_{K,\mathrm{design}}\,dx
=-\int_D S(u_1,\lambda_1)w\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

同时满足

\[
\int_D \theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0.
\]

矩阵形式为

\[
\begin{bmatrix}
A_\theta & c_\rho\\
c_\rho^T & 0
\end{bmatrix}
\begin{pmatrix}
\theta\\ p
\end{pmatrix}
=
\begin{pmatrix}
b_\rho\\0
\end{pmatrix}.
\]

其中

\[
(A_\theta\theta,w)
=
\int_D \nabla\theta:\nabla w\,dx
+
\int_D \theta\cdot w\,dx,
\]

\[
(c_\rho,w)=\int_D w\cdot\nabla\rho_{K,\mathrm{design}}\,dx,
\]

\[
(b_\rho,w)
=-\int_D S(u_1,\lambda_1)w\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

~~~freefem
Thetah [thetax, thetay], [Wx, Wy];

varf VelocityA([thetax, thetay], [Wx, Wy]) =
int2d(Th)(
        Grad(thetax, thetay)' * Grad(Wx, Wy)
      + [thetax, thetay]' * [Wx, Wy]
)
+ on(3, thetax = 0, thetay = 0)
+ on(11, 2, thetax = 0)
+ on(12, thetay = 0);

matrix Atheta = VelocityA(Thetah, Thetah);
set(Atheta, solver = sparsesolver);
cout << "Velocity matrix size = " << Atheta.n << " x " << Atheta.m << endl;

object << "Velocity matrix size = " << Atheta.n << " x " << Atheta.m << endl;


macro VelocityPDE(thetax, thetay, modex, modey, lambdaone, phiLS)
{
    varf Velocityc([thetax, thetay], [Wx, Wy]) =
    int2d(Th)(
        [Wx, Wy]' * gradRhoDesign(phiLS)
    )
    + on(3, thetax = 0, thetay = 0)
    + on(11, 2, thetax = 0)
    + on(12, thetay = 0);

    varf Velocityb([thetax, thetay], [Wx, Wy]) =
    - int2d(Th)(
        (
            2.0 * mu * (epsilon(modex, modey)' * epsilon(modex, modey))
          + lameLambda * div(modex, modey) * div(modex, modey)
          - lambdaone * (1.0 - weakM) / (1.0 - weak) * (modex^2 + modey^2)
        )
        * ([Wx, Wy]' * gradRhoDesign(phiLS))
    )
    + on(3, thetax = 0, thetay = 0)
    + on(11, 2, thetax = 0)
    + on(12, thetay = 0);

    real[int] cRho = Velocityc(0, Thetah);
    real[int] bRho = Velocityb(0, Thetah);

    real[int] zRho = Atheta^-1 * bRho;
    real[int] yRho = Atheta^-1 * cRho;

    real cy = cRho' * yRho;
    real cz = cRho' * zRho;
    real ptheta = 0.0;

    if (abs(cy) > 1e-30)
    {
        ptheta = cz / cy;
    }

    real[int] thetaSol = zRho - ptheta * yRho;

    thetax[] = thetaSol;
    cout << "Velocity Success! "<< endl;
} //
~~~

## 10. 监测量

程序中主要记录三个量：

1. 第一特征值

\[
\lambda_1;
\]

2. 当前体积

\[
\operatorname{Vol}(\phi)=\int_D\rho\,dx;
\]

3. 速度场对应的体积变化率

\[
d\operatorname{Vol}(\theta)=\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

如果速度场方程求解正确，则 `dVolTheta` 应该接近 \(0\)。

同时还记录

\[
d\lambda_1(\theta)
=-\int_D S(u_1,\lambda_1)\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

若该值大于 \(0\)，则说明当前速度场方向会使第一特征值上升。

~~~freefem
macro GetVol(Rho)
(
    int2d(Th)(Rho)
) //

macro GetdVol(thetax, thetay, phiLS)
(
    int2d(Th)([thetax, thetay]' * gradRhoDesign(phiLS))
) //

macro GetdLambda(modeX, modeY, thetax, thetay, lambda1, phiLS)
(
    - int2d(Th)(
        (
            2.0 * mu * (epsilon(modeX, modeY)' * epsilon(modeX, modeY))
          + lameLambda * div(modeX, modeY) * div(modeX, modeY)
          - lambda1 * (1.0 - weakM) / (1.0 - weak) * (modeX^2 + modeY^2)
        )
        * ([thetax, thetay]' * gradRhoDesign(phiLS))
      )
) //
~~~

## 11. dual-P^1 DG 水平集输运方程

水平集函数通过速度场 \(\theta\) 更新：

\[
\phi_t+\theta\cdot\nabla\phi=0.
\]

这里 \(\phi\) 使用间断有限元空间 `Phih=P1dc`。更新时使用隐式迎风 dual-P^1 DG 格式。若记 \(\phi^{n+1}\) 为下一步水平集函数，则单步离散大致对应

\[
\int_D \frac{\phi^{n+1}-\phi^n}{\Delta t}v\,dx
+
\int_D v\theta\cdot\nabla\phi^{n+1}\,dx
+
\text{upwind flux}=0.
\]

由于速度场在外边界上做了相应约束，可以减少界面从固定设计域边界非物理流出的情况。

~~~freefem
macro nn()
(
    N.x * thetax + N.y * thetay
) //

macro LFDG(phiOutput, phiInput, thetax, thetay, Tsum)
{
    real t = 0.0;
    real dt = Tsum / 5.0;

    Phih u, v;

    varf Adual(u, v) =
    int2d(Th)((u / dt + ([thetax, thetay]' * grad(u))) * v)
    - intalledges(Th)(
          (nTonEdge - 1.0) * v
        * (0.5 * abs(nn) - nn / 2.0) * jump(u)
      );

    varf DGRHS(u, v) = int2d(Th)((u / dt) * v);

    matrix Am = Adual(Phih, Phih);
    matrix Bm = DGRHS(Phih, Phih);
    set(Am, solver = sparsesolver);

    for (t = 0.0; t < Tsum; t += dt)
    {
        v[] = Bm * phiInput[];
        phiOutput[] = Am^-1 * v[];
        phiInput = phiOutput;
    }

    cout << "DG Success! "<< endl;
} //
~~~

## 12. 重新初始化

重新初始化使用 FreeFEM 的 `distance`。这里需要注意：

- 输入 `phiInput` 可以是一般有限元函数，例如当前的 `Phih` 中的 dual-P^1 DG 函数；
- 输出数组必须属于 `P1h`；
- 因此代码里先定义 `P1h phiRe`，然后执行

\[
\texttt{distance(Th, phiInput, phiRe[])}.
\]

最后再赋值回原来的水平集函数空间。

~~~freefem
macro RedistanceConvert(phiInput)
{
    P1h phiRe;
    distance(Th, phiInput, phiRe[]);
    phiInput = phiRe;
} //
~~~

## 13. 输出当前迭代信息

输出内容包括：

- 当前第一特征值 \(\lambda_1\)；
- 当前体积 \(\int_D\rho\,dx\)；
- 相比初始体积的体积变化；
- 当前速度场对应的体积变化率；
- 当前速度场对应的第一特征值变化率。

~~~freefem
macro OutPutData(i, eigValue, vol, dvol, dlambda, obj)
{
    if (i == 0)
    {
        obj << "Initial condition" << endl;
        obj << "lambda1 = " << eigValue[i] << endl;
        obj << "lambda2 = " << lambda2 << endl;
        obj << "gap12/lambda1 = " << gap12Ratio << endl;
        obj << "Vol = " << vol[i] << endl;
    }
    else
    {
        obj << i << "-th iteration" << endl;
        obj << "lambda1 = " << eigValue[i]
            << ";   diff = " << eigValue[i] - eigValue[i-1] << endl;
        obj << "Vol = " << vol[i]
            << ";   diff = " << vol[i] - vol[0] << endl;
        obj << "dVol(theta) = " << dvol << endl;
        obj << "dLambda(theta) = " << dlambda << endl;
        obj << "lambda2 = " << lambda2 << endl;
        obj << "gap12/lambda1 = " << gap12Ratio << endl;
    }
    obj.flush;
} //
~~~

## 14. 主迭代过程

每一次优化迭代执行以下步骤：

1. 已知当前 \(\phi^n\)，计算密度

\[
\rho_K^n=\chi_{D_{\mathrm{block}}}+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}(\phi^n),\qquad
\rho_M^n=\mathrm{blockMass}\,\chi_{D_{\mathrm{block}}}+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weakM}}(\phi^n).
\]

2. 解弹性结构第一特征值问题

\[
K(\rho_K^n)u_1^n=\lambda_1^nM(\rho_M^n)u_1^n.
\]

3. 由形状导数密度

\[
S(u_1^n,\lambda_1^n)
=
2\mu\varepsilon(u_1^n):\varepsilon(u_1^n)
+
\lambda_L(\operatorname{div}u_1^n)^2
-
\lambda_1^n\frac{1-\mathrm{weakM}}{1-\mathrm{weak}}|u_1^n|^2
\]

求保体积速度场 \(\theta^n\)。

4. 用 dual-P^1 DG 格式求解水平集输运方程

\[
\phi_t+\theta^n\cdot\nabla\phi=0.
\]

5. 每隔若干步做一次 `distance` 重新初始化，使 \(\phi\) 保持接近有向距离函数。

~~~freefem
Phih phitest;
real[int] EigenValue1(iterMax + 1);
real[int] EigenValue2(iterMax + 1);
real[int] EigenValue3(iterMax + 1);
real[int] EigenValue4(iterMax + 1);
real[int] Gap12Ratio(iterMax + 1);
real[int] Vol(iterMax + 1);
real dVolTheta, dLambdaTheta;

int iter = 0;

SolveFirstEigen(modeX, modeY, lambda1);
VelocityPDE(thetax, thetay, modeX, modeY, lambda1, phi);

EigenValue1[iter] = lambda1;
EigenValue2[iter] = lambda2;
EigenValue3[iter] = lambda3;
EigenValue4[iter] = lambda4;
Gap12Ratio[iter] = gap12Ratio;
Vol[iter] = GetVol(Rho);
dVolTheta = GetdVol(thetax, thetay, phi);
dLambdaTheta = GetdLambda(modeX, modeY, thetax, thetay, lambda1, phi);
OutPutData(iter, EigenValue1, Vol, dVolTheta, dLambdaTheta, object);

real tstart;
real tend;
real tcost;
real tsum = 0.0;
while (iter < iterMax)
{
    tstart = clock();
    cout << endl;
    cout << "Iter = " << iter << ", lambda1 = " << lambda1 << endl;
    object << endl << "Iter = " << iter << endl;

    phitest = phi;
    LFDG(phi, phitest, thetax, thetay, Titer);

    if ((iter - 1) % 8 == 0)
    {
        RedistanceConvert(phi);
        cout << "Redistance End!" << endl;
    }

    Rho = RhoWithBlock(phi);
    RhoMass = RhoMassWithBlock(phi);

    SolveFirstEigen(modeX, modeY, lambda1);
    VelocityPDE(thetax, thetay, modeX, modeY, lambda1, phi);

    iter++;
    EigenValue1[iter] = lambda1;
    EigenValue2[iter] = lambda2;
    EigenValue3[iter] = lambda3;
    EigenValue4[iter] = lambda4;
    Gap12Ratio[iter] = gap12Ratio;
    Vol[iter] = GetVol(Rho);
    dVolTheta = GetdVol(thetax, thetay, phi);
    dLambdaTheta = GetdLambda(modeX, modeY, thetax, thetay, lambda1, phi);

    object << "Titer = " << Titer << endl;
    OutPutData(iter, EigenValue1, Vol, dVolTheta, dLambdaTheta, object);

    plot(Rho, fill = 1, value = 1, cmm = "Rho" + iter, wait = false);
    if (iter % 10 == 0)
    {
        string filename = "./figure/Rho" + iter + ".vtu";
        savevtk(filename, Th, Rho);
        string filenameM = "./figure/RhoMass" + iter + ".vtu";
        savevtk(filenameM, Th, RhoMass);
    }

    tend = clock();
    tcost = tend - tstart;
    tsum = tsum + tcost;
    cout << "Iteration time = " << tcost << " s" << endl;
    object << "Iteration time = " << tcost << " s" << endl;
    object.flush;
}
cout << "Iteration sum time = " << tsum << " s" << endl;
object << "Iteration sum time = " << tsum << " s" << endl;
object.flush;
~~~

## 15. Matlab 数据输出

最终输出 Matlab 数据，包含每一步的体积、前四个特征值以及第一、第二特征值间隔比。

其中

\[
\texttt{Volume(i)}
\]

对应第 \(i\) 步的体积，

\[
\texttt{lambda1(i)}
\]

对应第 \(i\) 步的第一特征值；`gap12ratio(i)` 用于检测模态切换或近重特征值。

~~~freefem
ofstream OutMatlab(ModelMatlabName);
OutMatlab << "Num = " << iter + 1 << endl;
for (int i = 0; i <= iter; i++)
{
    OutMatlab << "Volume(" << i+1 << ") = " << Vol[i] << "; ";
    OutMatlab << "lambda1(" << i+1 << ") = " << EigenValue1[i] << "; ";
    OutMatlab << "lambda2(" << i+1 << ") = " << EigenValue2[i] << "; ";
    OutMatlab << "lambda3(" << i+1 << ") = " << EigenValue3[i] << "; ";
    OutMatlab << "lambda4(" << i+1 << ") = " << EigenValue4[i] << "; ";
    OutMatlab << "gap12ratio(" << i+1 << ") = " << Gap12Ratio[i] << ";" << endl;
}
~~~

## 16. 模型总结

这个程序实际求解的是如下优化问题：

\[
\max_{\phi}\lambda_1(\phi),
\]

其中 \(\lambda_1(\phi)\) 由广义特征值问题给出：

\[
\int_D
\rho_K(\phi)
\left(
2\mu\varepsilon(u):\varepsilon(v)
+
\lambda_L\operatorname{div}u\operatorname{div}v
\right)dx
=
\lambda
\int_D
\rho_M(\phi)u\cdot v\,dx.
\]

水平集函数通过

\[
\phi_t+\theta\cdot\nabla\phi=0
\]

演化。速度场 \(\theta\) 不是任意给定的，而是由保体积约束的下降方向方程计算：

\[
\int_D \nabla\theta:\nabla w\,dx
+
\int_D \theta\cdot w\,dx
+p\int_Dw\cdot\nabla\rho_{K,\mathrm{design}}\,dx
=
-\int_D S(u_1,\lambda_1)w\cdot\nabla\rho_{K,\mathrm{design}}\,dx,
\]

\[
\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0.
\]

其中

\[
S(u_1,\lambda_1)
=
2\mu\varepsilon(u_1):\varepsilon(u_1)
+
\lambda_L(\operatorname{div}u_1)^2
-
\lambda_1\frac{1-\mathrm{weakM}}{1-\mathrm{weak}}|u_1|^2.
\]

由于速度场满足

\[
\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0,
\]

所以体积在连续层面满足

\[
\frac{d}{dt}\int_D\rho\,dx=0.
\]

这就是该程序不需要额外添加任意体积乘子的原因。体积约束已经包含在速度场方程中。
