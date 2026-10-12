# 弹性结构第一特征值水平集优化

本文档把数学模型说明和 FreeFEM 代码融合在一起。所有 `~~~freefem` 代码块从上到下拼接后，就是一个完整的 `.edp` 程序。

本程序在原来的水平集特征值拓扑优化框架上，把计算模型改为论文 Figure 6 中的双端固支梁模型。设计域为 \(4 \times 1\) 的长方形按比例缩小为 \(3.2 \times 0.8\)，也就是缩放因子为 \(0.8\)。结构位移在左、右端满足 Dirichlet 边界条件，上、下边界为自然边界条件。速度场只在整个外边界上满足无穿透条件 \(\theta\cdot n=0\)。中心小黑块作为非设计附加质量处理，质量贡献约为全域实体质量的 \(25\%\)。

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
| `Lx`, `Ly` | 设计域长和高，这里为 `3.2` 和 `0.8` |
| `xmin`, `xmax`, `ymin`, `ymax` | 设计域边界坐标 |
| `blockHalf` | 中心小黑块的半边长 |
| `blockMassRatio` | 中心小黑块质量贡献相对于全域实体质量的比例，这里取 `0.25` |
| `blockMass` | 中心小黑块在质量矩阵中的等效密度 |
| `mu` | 线弹性 Lamé 常数 \(\mu\) |
| `lameLambda` | 线弹性 Lamé 常数 \(\lambda_L\)。这里不用 `lambda` 作为变量名，避免与特征值概念混淆 |
| `weak` | 刚度矩阵中的弱材料系数，即 \(\rho_K=\mathrm{weak}+(1-\mathrm{weak})H_\gamma(\phi)\) |
| `weakM` | 质量矩阵中的弱材料系数。这里取 `1e-4`，即空洞不提供质量 |
| `gamma` | 光滑 Heaviside 函数的界面宽度 |
| `sigma` | 特征值求解的 shift 参数。这里取 `0.0` |
| `nev` | 计算的特征值个数。虽然当前只优化第一个特征值，但多算几个有利于观察模态交换 |
| `object` | 优化过程信息输出文件 |
| `ModelMatlabName` | 输出给 Matlab 画图的数据文件 |

目标是最大化第一特征值 \(\lambda_1\)。实际速度场方程按最小化

\[
J(\phi)=-\lambda_1(\phi)
\]

来写，因此速度场方向是使 \(J\) 下降，也就是使 \(\lambda_1\) 上升的方向。

中心小黑块设为

\[
D_{\mathrm{block}}=[-0.08,0.08]\times[-0.08,0.08].
\]

中心小黑块取为正方形区域 \([-0.08,0.08]\times[-0.08,0.08]\)，边长为 \(0.16\)，约为当前梁高 \(0.6\) 的 \(25\%\)。

由于 \(|D|=3.2\times0.8=2.56\)，\(|D_{\mathrm{block}}|=0.16^2=0.0256\)，为了使中心小黑块质量约占全域实体质量的 \(25\%\)，取

\[
\mathrm{blockMass}=\frac{0.25|D|}{|D_{\mathrm{block}}|}=25.
\]

~~~freefem
real Titer = 0.03;
int iterMax = 1000;

real Lx = 3.2;
real Ly = 0.8;
real xmin = -Lx / 2.0;
real xmax =  Lx / 2.0;
real ymin = -Ly / 2.0;
real ymax =  Ly / 2.0;

real blockHalf = 0.08;
real blockMassRatio = 0.25;
real domainMass = Lx * Ly;
real blockArea = (2.0 * blockHalf) * (2.0 * blockHalf);
real blockMass = blockMassRatio * domainMass / blockArea;

real mu = 8.0;
real lameLambda = 1.0;
real weak = 1e-2;
real weakM = 1e-4;
real gamma = 0.005;

real sigma = 0.0;
int nev = 4;
real modeGapTol = 1e-2;

ofstream object("./Output/MainEigen.obj");
string ModelMatlabName = "matlab/Model1EigenShape2.m";

object << "Lx = " << Lx << ", Ly = " << Ly << endl;
object << "blockMass = " << blockMass << endl;
object << "block mass ratio = " << blockMass * blockArea / domainMass << endl;
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

本程序取

\[
\mathrm{weakM}=0,
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
    (x >= -blockHalf && x <= blockHalf && y >= -blockHalf && y <= blockHalf)
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

中心小黑块区域为

\[
D_{\mathrm{block}}=[-0.08,0.08]\times[-0.08,0.08].
\]


该区域作为非设计区域处理，但刚度和质量的含义不同。刚度矩阵中，小黑块仍是普通实心材料，

\[
\rho_K(x)=1,\qquad x\in D_{\mathrm{block}}.
\]

质量矩阵中，小黑块用于模拟中心附加质量，因此取

\[
\rho_M(x)=\mathrm{blockMass},\qquad x\in D_{\mathrm{block}}.
\]

在其余可设计区域中，刚度密度取 \(H_{\gamma,\mathrm{weak}}(\phi)\)，质量密度取 \(H_{\gamma,\mathrm{weakM}}(\phi)\)，且这里 \(\mathrm{weakM}=0\)。因此实际使用

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
D=(-2.1,2.1)\times(-0.3,0.3),
\]

也就是 Figure 6 中 \(14\times2\) 设计域按比例缩小后的 \(4.2\times0.6\) 长方形。

边界标号改为：

| label | 边界 | 弹性边界条件 | 速度场边界条件 |
|:---:|:---|:---|:---|
| `1` | 左边界 \(x=-2.1\) | \(u=0\) | \(\theta_x=0\) |
| `2` | 右边界 \(x=2.1\) | \(u=0\) | \(\theta_x=0\) |
| `3` | 上边界 \(y=0.3\) | 自然边界条件 | \(\theta_y=0\) |
| `4` | 下边界 \(y=-0.3\) | 自然边界条件 | \(\theta_y=0\) |

因此弹性结构只有左、右端为 Dirichlet 边界，上、下边界不加 `on()`，即自动对应自然边界条件。速度场在整个边界满足

\[
\theta\cdot n=0.
\]

有限元空间选取为：

| 空间 | 有限元 | 作用 |
|:---:|:---:|:---|
| `Uh` | \([P2,P2]\) | 弹性特征振型空间 |
| `Thetah` | \([P1,P1]\) | 水平集速度场空间 |
| `P1h` | \(P1\) | `distance` 重新初始化的输出空间 |
| `Phih` | \(P1dc\) | dual-P^1 DG 水平集函数空间 |

这里要注意：`distance` 可以接受一般有限元函数或者表达式作为输入，但输出数组必须属于 `P1h`。因此重新初始化时先输出到 `P1h phiRe`，再赋值回 `Phih`。

~~~freefem
border bdbottom(t = xmin, xmax) { x=t;    y=ymin; label=4; };
border bdright(t = ymin, ymax)  { x=xmax; y=t;    label=2; };
border bdup(t = xmax, xmin)     { x=t;    y=ymax; label=3; };
border bdleft(t = ymax, ymin)   { x=xmin; y=t;    label=1; };

int n = 1;
mesh Th = buildmesh(
      bdbottom(200*n) + bdright(50*n)
    + bdup(200*n) + bdleft(50*n)
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

初始材料区域参考论文 Figure 7(a) 的孔洞排布方式。水平集符号约定为

\[
\phi_0>0
\]

对应孔洞区域，\(\phi_0<0\) 对应材料区域。之后通过 `distance` 重新初始化为近似有向距离函数。由于密度函数中约定 \(\phi<0\) 为实心材料，所以这个符号约定与后续 `Heaviside` 一致。

圆孔半径取为两位小数

\[
r_h=0.06.
\]

从左到右共布置 13 列孔洞。第 1 列和第 11 列贴在左右边界上，分别由上角四分之一圆孔、中间贴边半圆孔和下角四分之一圆孔组成。去掉左右两端之后，中间 11 列按 2 孔列和 3 孔列交错排列：2 孔列包含上下两个内部圆孔；3 孔列包含上边界半圆孔、中间内部圆孔和下边界半圆孔。中心小黑块是非设计区域，后续密度函数会强制该区域为实心材料和附加质量区域。

~~~freefem
real rHole = 0.12;
real xStepHole = Lx / 14.0;
real yTwoHole = Ly / 4.0;

P1h phi0;
phi0 = -1.0;

for (int i = 0; i <= 18; i++)
{
    real xc = xmin + i * xStepHole;

    if (i == 0 || i == 18)
    {
        phi0 = max(phi0, rHole - sqrt((x - xc)^2 + (y - ymax)^2));
        phi0 = max(phi0, rHole - sqrt((x - xc)^2 + y^2));
        phi0 = max(phi0, rHole - sqrt((x - xc)^2 + (y - ymin)^2));
    }
    else
    {
        if (i % 2 == 1)
        {
            phi0 = max(phi0, rHole - sqrt((x - xc)^2 + (y - yTwoHole)^2));
            phi0 = max(phi0, rHole - sqrt((x - xc)^2 + (y + yTwoHole)^2));
        }
        else
        {
            phi0 = max(phi0, rHole - sqrt((x - xc)^2 + (y - ymax)^2));
            phi0 = max(phi0, rHole - sqrt((x - xc)^2 + y^2));
            phi0 = max(phi0, rHole - sqrt((x - xc)^2 + (y - ymin)^2));
        }
    }
}

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

这里取 `sigma=0`，所以它就是刚度矩阵 \(K(
\rho)\)。Dirichlet 边界条件只施加在左、右端，即 `label=1` 和 `label=2`。上、下边界不加 `on()`，因此对应自然边界条件。刚度矩阵和质量矩阵都使用同一组 Dirichlet 自由度约束。实际在 FreeFEM 中，为避免边界自由度产生靠近目标特征值的伪特征值，刚度型矩阵和质量矩阵使用不同的 `tgv`：刚度矩阵取很大的 `tgv=1e30`，质量矩阵取 `tgv=1e10`。这样边界自由度对应的伪特征值非常大，不会影响靠近零的第一正特征值。

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
+ on(1, 2, ux = 0, uy = 0);

varf Mass([ux, uy], [vx, vy]) =
int2d(Th)(RhoMass * (ux * vx + uy * vy))
+ on(1, 2, ux = 0, uy = 0);

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

这里按当前代码中的 ersatz 密度归一化方式，把质量项写为 \(\lambda_1(1-\mathrm{weakM})|u_1|^2/(1-\mathrm{weak})\)。由于本程序取 \(\mathrm{weakM}=0\)，代码中等价于 \(\lambda_1|u_1|^2/(1-\mathrm{weak})\)。这一点会同步体现在速度场右端和 `dLambdaTheta` 的监测量中。

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

这里采用速度场强制体积守恒，所以不需要另外在目标函数中手动加入任意体积乘子。体积为

\[
\operatorname{Vol}(\phi)=\int_D \rho(\phi)\,dx.
\]

沿水平集输运方程演化时，小黑块部分固定不变，所以

\[
\frac{d}{dt}\operatorname{Vol}(\phi)=\int_D \rho_t\,dx=-\int_D \theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

因此只要速度场满足

\[
\int_D \theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0,
\]

就得到连续层面的体积守恒。

本模型中速度场边界条件只取无穿透条件：

\[
\theta\cdot n=0\qquad\text{on }\partial D.
\]

对当前矩形区域，这等价于

\[
\theta_x=0\quad\text{on }x=\pm2.1,
\qquad
\theta_y=0\quad\text{on }y=\pm0.3.
\]

于是求速度场 \(\theta\) 和约束变量 \(p\)，使得

\[
\int_D \nabla\theta:\nabla w\,dx+
\int_D \theta\cdot w\,dx
+p\int_D w\cdot\nabla\rho_{K,\mathrm{design}}\,dx
=-dJ(w),
\]

即

\[
\int_D \nabla\theta:\nabla w\,dx+
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
(A_\theta\theta,w)=\int_D \nabla\theta:\nabla w\,dx+\int_D \theta\cdot w\,dx,
\]

\[
(c_\rho,w)=\int_D w\cdot\nabla\rho_{K,\mathrm{design}}\,dx,
\]

\[
(b_\rho,w)=-\int_D S(u_1,\lambda_1)w\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

~~~freefem
Thetah [thetax, thetay], [Wx, Wy];

varf VelocityA([thetax, thetay], [Wx, Wy]) =
int2d(Th)(
        Grad(thetax, thetay)' * Grad(Wx, Wy)
      + [thetax, thetay]' * [Wx, Wy]
)
+ on(1, 2, thetax = 0)
+ on(3, 4, thetay = 0);

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
    + on(1, 2, thetax = 0)
    + on(3, 4, thetay = 0);

    varf Velocityb([thetax, thetay], [Wx, Wy]) =
    - int2d(Th)(
        (
            2.0 * mu * (epsilon(modex, modey)' * epsilon(modex, modey))
          + lameLambda * div(modex, modey) * div(modex, modey)
          - lambdaone * (1.0 - weakM) / (1.0 - weak) * (modex^2 + modey^2)
        )
        * ([Wx, Wy]' * gradRhoDesign(phiLS))
    )
    + on(1, 2, thetax = 0)
    + on(3, 4, thetay = 0);

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

由于速度场在整个外边界满足 \(\theta\cdot n=0\)，水平集输运在连续层面不会通过外边界产生通量。

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
    if (iter <= 150)
        Titer = 0.03;
    else if (150 < iter && iter <= 300)
        Titer = 0.02;
    else if (300 < iter && iter <= 450)
        Titer = 0.01;
    else if (450 < iter && iter <= 600)
        Titer = 0.008;
    else 
        Titer = 0.004;
    tstart = clock();
    cout << endl;
    cout << "Iter = " << iter << ", lambda1 = " << lambda1 << endl;
    object << endl << "Iter = " << iter << endl;

    phitest = phi;
    LFDG(phi, phitest, thetax, thetay, Titer);

    if ((iter - 1) % 20 == 0)
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

    plot(Rho, fill = 1, value = 0, cmm = "Rho" + iter, wait = false);
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

其中设计域为

\[
D=(-2.1,2.1)\times(-0.3,0.3).
\]

弹性边界条件为

\[
u=0\quad\text{on }\Gamma_D=\{x=-2.1\}\cup\{x=2.1\},
\]

上、下边界为自然边界条件：

\[
\sigma(u)n=0\quad\text{on }\Gamma_N=\{y=-0.3\}\cup\{y=0.3\}.
\]

第一特征值由广义特征值问题给出：

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

中心小黑块为

\[
D_{\mathrm{block}}=[-0.05,0.05]\times[-0.05,0.05],
\]

它不参与优化。刚度矩阵中小黑块为普通实心材料，质量矩阵中小黑块具有等效密度 `blockMass = 40.32`，从而

\[
\mathrm{blockMass}\,|D_{\mathrm{block}}|
\approx 0.16|D|.
\]

水平集函数通过

\[
\phi_t+\theta\cdot\nabla\phi=0
\]

演化。速度场在整个外边界满足无穿透条件

\[
\theta\cdot n=0\quad\text{on }\partial D,
\]

也就是左、右边界取 \(\theta_x=0\)，上、下边界取 \(\theta_y=0\)。速度场由保体积约束的下降方向方程计算：

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

