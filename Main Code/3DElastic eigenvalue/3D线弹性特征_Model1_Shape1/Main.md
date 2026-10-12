# 三维弹性结构第一特征值水平集优化

本文档把数学模型说明和 FreeFEM 代码融合在一起。所有 `~~~freefem` 代码块从上到下拼接后，就是一个完整的 `.edp` 程序。

本程序考虑三维悬臂梁的第一特征值最大化问题。模型设置参考文献中 Figure 7.7 的三维算例：设计区域为

\[
D=(0,1)\times(0,2)\times(0,1),
\]

左侧面施加零位移边界条件，即在与 \(D_{\mathrm{block}}\) 正对的 \(x=0\) 的 \(yz\) 平面上固定。右侧面中心有一个小质量块区域，该区域不参与优化，且质量密度取为普通材料的 `333.0` 倍，使该块质量约占整体质量的 25%。初始构型采用带周期孔的三维悬臂梁，孔沿厚度方向贯穿，用来近似 Figure 7.7 中的三维初始化。

优化目标是最大化第一特征值

\[
\lambda_1=\omega_1^2.
\]

实际速度场按最小化

\[
J(\phi)=-\lambda_1(\phi)
\]

来构造。速度场通过带体积约束的 Riesz 表示求得，所以不需要在目标函数中额外手动加入任意体积乘子。

## 1. 必要的 load

三维模型需要使用三维网格模块、VTK 输出模块、水平集重新初始化模块以及 PETSc 模块。

- `msh3` 用于三维网格；
- `tetgen` 用于由边界曲面生成三维四面体网格；
- `distance` 用于把水平集函数重新初始化为近似有向距离函数；
- `iovtk` 用于输出 `.vtu` 文件，方便 ParaView 后处理；
- `PETSc` 用于 `Mat` 类型的线性系统写法，本文档中三维 dual-P^1 DG 水平集输运方程采用 PETSc 的 `Mat` 矩阵表示；
- `MUMPS_seq` 提供串行 MUMPS 稀疏直接求解器；
- 弹性特征值问题和速度场线性系统都把矩阵求解器设为 `sparsesolver`。加载 `MUMPS_seq` 后，`sparsesolver` 使用串行 MUMPS 后端，避免使用 MPI 并行。

~~~freefem
load "msh3"
load "tetgen"
load "distance"
load "iovtk"
load "PETSc"
load "MUMPS_seq"
verbosity = 0;
~~~

## 2. 参数设置

| 变量名 | 含义 |
|:---:|:---|
| `Titer` | 每次形状更新的伪时间总步长 |
| `iterMax` | 最大优化迭代步数 |
| `hMeshTarget` | 三维网格目标尺寸，代码中取 `0.09` |
| `Nx, Ny, Nz` | 保留的结构网格剖分参数，当前 tetgen 网格生成未直接使用 |
| `Young` | 杨氏模量，按照文献算例归一化取 `1.0` |
| `nu` | Poisson 比，按照文献算例取 `0.3` |
| `mu` | Lamé 常数 \(\mu=E/(2(1+\nu))\) |
| `lameLambda` | Lamé 常数 \(\lambda_L=E\nu/((1+\nu)(1-2\nu))\) |
| `weak` | 刚度矩阵中的弱材料系数，取 `1e-2` |
| `weakM` | 质量矩阵中的弱材料系数，取 `1e-4` |
| `blockMass` | 右端小质量块的质量密度，取普通材料的 `333.0` 倍 |
| `gamma` | 光滑 Heaviside 函数的界面宽度，这里取为一个网格尺寸量级 |
| `sigma` | 特征值求解 shift 参数，用于寻找靠近 `0` 的最小正特征值 |
| `nev` | 计算的特征值个数，用于检测模态切换；当前取 `6` |
| `modeGapTol` | 第一、第二特征值相对间隔过小时的提示阈值 |
| `object` | 优化过程信息输出文件 |
| `ModelMatlabName` | 输出给 Matlab 画图的数据文件 |

文献中的三维算例使用 \(1\times2\times1\) 的工作区域，并在 \(15\times30\times15\) 网格上计算。这里改用 `tetgen` 生成四面体网格：代码先取

\[
\texttt{meshEdgeTarget}=0.64\,\texttt{hMeshTarget},
\]

再由该目标边长计算外边界面的剖分数 `surfNx`、`surfNy`、`surfNz`，最后用 `tetg` 生成三维网格。

~~~freefem
real Titer = 0.02;
int iterMax = 300;

real hMeshTarget = 0.09;
int Nx = 20;
int Ny = 40;
int Nz = 20;

real Young = 1.0;
real nu = 0.3;
real mu = Young / (2.0 * (1.0 + nu));
real lameLambda = Young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

real weak = 1e-2;
real weakM = 1e-4;
real blockMass = 333.0;
real gamma = hMeshTarget;

real sigma = 0.0;
int nev = 6;
real modeGapTol = 1e-2;

ofstream object("./Output/MainEigen3D.obj");
string ModelMatlabName = "matlab/ModelEigen3D.m";
~~~

## 3. 水平集函数和密度函数

材料区域由水平集函数表示：

\[
\Omega=\{x\in D:\phi(x)<0\}.
\]

空洞或弱材料区域对应

\[
\phi(x)>0.
\]

为了在固定网格上求解弹性广义特征值问题，引入刚度密度和质量密度。代码中的 Heaviside 函数直接返回密度插值。刚度矩阵使用

\[
\rho_K(\phi)=H_{\gamma,\mathrm{weak}}(\phi),
\]

质量矩阵使用

\[
\rho_M(\phi)=H_{\gamma,\mathrm{weakM}}(\phi).
\]

这里的 \(H_\gamma\) 按照 \(\phi<0\) 为材料、\(\phi>0\) 为空洞的约定写成：

\[
H_{\gamma,a}(\phi)
=
\begin{cases}
a, & \phi>\gamma,\\[1mm]
(1-a)\left(\dfrac12-\dfrac{\phi}{\gamma}+\dfrac{\phi^2}{2\gamma^2}\right)+a,
&0<\phi\le \gamma,\\[2mm]
(1-a)\left(\dfrac12-\dfrac{\phi}{\gamma}-\dfrac{\phi^2}{2\gamma^2}\right)+a,
&-\gamma\le \phi\le 0,\\[2mm]
1, & \phi< -\gamma.
\end{cases}
\]

其中刚度密度取 \(a=\mathrm{weak}\)，质量密度取 \(a=\mathrm{weakM}\)。

右端中心小质量块区域记为

\[
D_{\mathrm{block}}
=
\{(x,y,z):0.9\le x\le 1,
\ 0.9\le y\le 1.1,
\ 0.45\le z\le 0.55\}.
\]

该块尺寸为 \(0.1\times0.2\times0.1\)，对应 Figure 7.7 中右端中心的非设计重质量区域。

该区域不参与优化。刚度矩阵中它保持为普通实心材料，质量矩阵中它的密度为普通材料的 `333.0` 倍：

\[
\rho_K(x)=1,\qquad x\in D_{\mathrm{block}},
\]

\[
\rho_M(x)=333.0,\qquad x\in D_{\mathrm{block}}.
\]

由于 \(|D_{\mathrm{block}}|=0.1\times0.2\times0.1=0.002\)，而 \(|D|=2\)，若其余区域按单位密度估算，则
\[
\frac{333.0\times0.002}{333.0\times0.002+(2-0.002)}
\approx 25\%.
\]
因此 `blockMass = 333.0` 对应小质量块质量约占整体质量的 25%。

因此代码中实际使用

\[
\rho_K(\phi)
=
\chi_{D_{\mathrm{block}}}
+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}(\phi),
\]

\[
\rho_M(\phi)
=
\mathrm{blockMass}\,\chi_{D_{\mathrm{block}}}
+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weakM}}(\phi).
\]

形状导数和体积约束只作用在可设计区域，所以使用

\[
\nabla\rho_{K,\mathrm{design}}
=(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}'(\phi)\nabla\phi.
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

macro Hprime(phi)
(
    (0 < phi && phi <= gamma)
    * ((1.0 - weak) * (-1.0/gamma + phi/(gamma^2)))
  + (-gamma <= phi && phi <= 0)
    * ((1.0 - weak) * (-1.0/gamma - phi/(gamma^2)))
) //

macro BlackBlock()
(
    (x >= 0.9 && x <= 1.0
  && y >= 0.9 && y <= 1.1
  && z >= 0.45 && z <= 0.55)
) //

macro RhoWithBlock(phi)
(
    BlackBlock * 1.0
  + (1.0 - BlackBlock) * HeaviSide(phi)
) //

macro RhoMassWithBlock(phi)
(
    BlackBlock * blockMass
  + (1.0 - BlackBlock)
    * (
          (phi > gamma) * weakM
        + (0 < phi && phi <= gamma)
          * ((1.0 - weakM) * (0.5 - phi/gamma + 0.5*phi^2/gamma^2) + weakM)
        + (-gamma <= phi && phi <= 0)
          * ((1.0 - weakM) * (0.5 - phi/gamma - 0.5*phi^2/gamma^2) + weakM)
        + (phi < -gamma) * 1.0
      )
) //
~~~

## 4. 三维线弹性记号

位移场为

\[
u=(u_x,u_y,u_z).
\]

应变张量为

\[
\varepsilon(u)=\frac12(\nabla u+\nabla u^T).
\]

在三维情形下，代码中使用带 \(\sqrt2\) 权重的 Voigt 形式：

\[
\varepsilon(u)
=
\left(
\partial_xu_x,
\partial_yu_y,
\partial_zu_z,
\frac{\partial_yu_x+\partial_xu_y}{\sqrt2},
\frac{\partial_zu_x+\partial_xu_z}{\sqrt2},
\frac{\partial_zu_y+\partial_yu_z}{\sqrt2}
\right).
\]

这样代码中的向量内积正好对应张量双点积

\[
\varepsilon(u):\varepsilon(v).
\]

散度为

\[
\operatorname{div}u=\partial_xu_x+\partial_yu_y+\partial_zu_z.
\]

~~~freefem
macro grad(u) [dx(u), dy(u), dz(u)] //

macro epsilon(ux, uy, uz)
[
    dx(ux),
    dy(uy),
    dz(uz),
    (dy(ux) + dx(uy))/sqrt(2.0),
    (dz(ux) + dx(uz))/sqrt(2.0),
    (dz(uy) + dy(uz))/sqrt(2.0)
] //

macro Grad(ux, uy, uz)
[
    dx(ux), dy(ux), dz(ux),
    dx(uy), dy(uy), dz(uy),
    dx(uz), dy(uz), dz(uz)
] //

macro div(ux, uy, uz) (dx(ux) + dy(uy) + dz(uz)) //

macro gradRhoDesign(phi)
(
    (1.0 - BlackBlock) * Hprime(phi) * grad(phi)
) //
~~~

## 5. 三维网格和有限元空间

计算区域为

\[
D=(0,1)\times(0,2)\times(0,1).
\]

边界标号设置为：

| label | 边界位置 | 用途 |
|:---:|:---|:---|
| `2` | \(x=0\) 的 \(yz\) 平面 | 与 \(D_{\mathrm{block}}\) 正对的左侧固定端，弹性位移满足 \(u=0\)，速度场满足 \(\theta_x=0\) |
| `3` | \(x=1\) 右端中心小块对应的边界面 | 与小质量块相交的边界，速度场三方向取零 |
| `11` | \(x=1\) 的其余部分 | 右端外边界，约束 \(\theta_x=0\) |
| `12` | \(y=0\) 或 \(y=2\) | 上下外边界，约束 \(\theta_y=0\) |
| `13` | \(z=0\) 或 \(z=1\) | 前后外边界，约束 \(\theta_z=0\) |

因此，若边界恰好落在 block 区域上，速度场使用三方向齐次 Dirichlet 条件；其余外边界只约束法向速度为零。左端虽然是弹性固定端，但速度场仍按滑移边界处理，即在 \(x=0\) 上只约束 \(\theta\cdot n=0\)，对应 \(\theta_x=0\)。

三维特征值问题的自由度比二维大很多，因此这里使用线性有限元近似弹性振型。这与文献中使用低阶单元进行三维计算的思路一致，也能避免网格较细时系统维数过大。

| 空间 | 有限元 | 作用 |
|:---:|:---:|:---|
| `Uh` | \([P1,P1,P1]\) | 三维弹性特征振型空间 |
| `Thetah` | \([P1,P1,P1]\) | 三维水平集速度场空间 |
| `P1h` | \(P1\) | `distance` 重新初始化的输出空间 |
| `Phih` | \(P1dc\) | dual-P^1 DG 水平集函数空间 |

`distance` 的输出数组必须属于连续的 `P1h` 空间。因此每次重新初始化时先得到 `P1h phiRe`，再赋值回 `Phih`。

~~~freefem
real epsLabel = 1e-8;

func Labelfunc =
      (x < epsLabel) ? 2
    : ((x > 1.0 - epsLabel
        && y >= 0.9 && y <= 1.1
        && z >= 0.45 && z <= 0.55) ? 3
    : ((x > 1.0 - epsLabel) ? 11
    : (((y < epsLabel) || (y > 2.0 - epsLabel)) ? 12
    : (((z < epsLabel) || (z > 1.0 - epsLabel)) ? 13 : 14))));

real meshEdgeTarget = 0.64 * hMeshTarget;
int surfNx = max(2, int(ceil(1.0 / meshEdgeTarget)));
int surfNy = max(2, int(ceil(2.0 / meshEdgeTarget)));
int surfNz = max(2, int(ceil(1.0 / meshEdgeTarget)));

mesh FaceXY = square(surfNx, surfNy, [x, 2.0*y]);
mesh FaceXZ = square(surfNx, surfNz, [x, y]);
mesh FaceYZ = square(surfNy, surfNz, [2.0*x, y]);

meshS SurfZ0 = movemesh23(FaceXY, transfo = [x, y, 0.0], orientation = -1);
meshS SurfZ1 = movemesh23(FaceXY, transfo = [x, y, 1.0], orientation = 1);
meshS SurfY0 = movemesh23(FaceXZ, transfo = [x, 0.0, y], orientation = 1);
meshS SurfY1 = movemesh23(FaceXZ, transfo = [x, 2.0, y], orientation = -1);
meshS SurfX0 = movemesh23(FaceYZ, transfo = [0.0, x, y], orientation = -1);
meshS SurfX1 = movemesh23(FaceYZ, transfo = [1.0, x, y], orientation = 1);

meshS BoxSurf = SurfZ0 + SurfZ1 + SurfY0 + SurfY1 + SurfX0 + SurfX1;
BoxSurf = change(BoxSurf, rmInternalEdges = true);

real tetVolume = meshEdgeTarget^3 / 6.0;
real[int] tetRegion = [0.5, 1.0, 0.5, 1, tetVolume];
mesh3 ThRaw = tetg(BoxSurf, switch = "pqaAAYYQ", nbofregions = 1, regionlist = tetRegion);
mesh3 Th = change(ThRaw, flabel = Labelfunc);

fespace Uh(Th, [P1, P1, P1]);
fespace Thetah(Th, [P1, P1, P1]);
fespace P1h(Th, P1);
fespace Phih(Th, P1dc);

P1h h1 = hTriangle;
object << "maximal size of a mesh tetrahedron h = " << h1[].max << endl;
object << "minimal size of a mesh tetrahedron h = " << h1[].min << endl;
cout << "maximal size of a mesh tetrahedron h = " << h1[].max << endl;
cout << "minimal size of a mesh tetrahedron h = " << h1[].min << endl;
~~~

## 6. 初始形状

初始材料区域参考 Figure 7.7 的三维悬臂梁初始化：整体为一根三维悬臂梁，内部布置周期孔。这里把初始孔近似为沿 \(z\) 方向贯穿的椭圆柱孔。令

\[
\phi_0(x,y,z)>0
\]

表示孔洞区域，

\[
\phi_0(x,y,z)<0
\]

表示材料区域。

为了得到多个孔洞，先在 `P1h` 空间中构造所有椭圆柱孔的并集：

\[
\phi_{\mathrm{hole}}(x,y,z)
=
\max_{i,j}\left\{
1-
\sqrt{\left(\frac{x-x_i}{\mathrm{holeRx}}\right)^2
+\left(\frac{y-y_j}{\mathrm{holeRy}}\right)^2}
\right\}.
\]

这里孔沿 \(z\) 方向贯穿，所以表达式不含 \(z\)。这会在三维长方体的前后面形成周期孔结构，并在侧面形成截断孔结构。代码没有在孔阵列生成阶段专门跳过 block 区域；而是在密度函数 `RhoWithBlock` 和 `RhoMassWithBlock` 中通过 `BlackBlock` 强制 \(D_{\mathrm{block}}\) 保持普通实心刚度和重质量密度。

然后使用 `distance` 重新初始化：

\[
\texttt{distance(Th, phiHole, phiInit[])}.
\]

最后对 `phiInit` 做常数平移校准，使初始刚度密度体积接近整体区域体积的 50%，然后计算初始刚度密度和质量密度。

校准形式为

\[
\phi=\phi_{\mathrm{init}}+s,
\]

其中 \(s\) 由二分搜索确定，使 \(\int_D \rho_K(\phi)\,dx\) 接近 `InitialVolumeRatio * int3d(Th)(1.0)`。默认 `InitialVolumeRatio = 0.5`。

~~~freefem
P1h phiHole = -1.0;
real holeRx = 0.095;
real holeRy = 0.145;

for (int j = 0; j < 6; j++)
{
    real holeYC = 0.4*j;

    phiHole = max(phiHole,
                  1.0 - sqrt(((x - 0.0)/holeRx)^2
                             + ((y - holeYC)/holeRy)^2));

    phiHole = max(phiHole,
                  1.0 - sqrt(((x - 1.0)/holeRx)^2
                             + ((y - holeYC)/holeRy)^2));
}

for (int j = 0; j < 5; j++)
{
    real holeYC = 0.2 + 0.4*j;

    phiHole = max(phiHole,
                  1.0 - sqrt(((x - 0.5)/holeRx)^2
                             + ((y - holeYC)/holeRy)^2));
}

P1h phiInit;
distance(Th, phiHole, phiInit[]);

Phih phi, Rho, RhoMass;
phi = phiInit;

real domainVolume = int3d(Th)(1.0);
real targetInitialVolume = initialVolumeRatio * domainVolume;
real initialPhiShiftLow = -2.0;
real initialPhiShiftHigh = 2.0;
real initialPhiShift = 0.0;
real initialCalibratedVolume = 0.0;

for (int k = 0; k < 50; k++)
{
    initialPhiShift = 0.5 * (initialPhiShiftLow + initialPhiShiftHigh);
    phi = phiInit + initialPhiShift;
    Rho = RhoWithBlock(phi);
    initialCalibratedVolume = int3d(Th)(Rho);

    if (initialCalibratedVolume > targetInitialVolume)
    {
        initialPhiShiftLow = initialPhiShift;
    }
    else
    {
        initialPhiShiftHigh = initialPhiShift;
    }
}

initialPhiShift = 0.5 * (initialPhiShiftLow + initialPhiShiftHigh);
phi = phiInit + initialPhiShift;
Rho = RhoWithBlock(phi);
RhoMass = RhoMassWithBlock(phi);

savevtk("./figure/Rho0_3D.vtu", Th, Rho);
savevtk("./figure/phi0_3D.vtu", Th, phi);
~~~

## 7. 弹性结构第一特征值问题

给定当前水平集函数 \(\phi\)，求三维广义特征值问题

\[
K(\rho_K)u=\lambda M(\rho_M)u.
\]

弱形式为：求 \((\lambda,u)\)，使得

\[
a_{\rho_K}(u,v)=\lambda m_{\rho_M}(u,v),
\qquad \forall v\in [H^1_{\Gamma_D}(D)]^3.
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

左侧固定端，即与 \(D_{\mathrm{block}}\) 正对的 \(x=0\) 的 \(yz\) 平面，满足

\[
u=0\qquad \text{on }\{x=0\}.
\]

在 FreeFEM 中，为避免固定边界自由度产生靠近目标特征值的伪特征值，刚度矩阵使用很大的 `tgv`，质量矩阵使用普通 `tgv`。这样固定边界自由度对应的伪特征值会远离靠近零的第一正特征值。

~~~freefem
Uh [modeX, modeY, modeZ], [testX, testY, testZ];
real lambda1, lambda2, lambda3, lambda4;
real lambda5, lambda6;
real gap12, gap12Ratio;
real lastEigenSolveTime = 0.0;
real lastVelocitySolveTime = 0.0;

varf EigenOP([ux, uy, uz], [vx, vy, vz]) =
int3d(Th)(
      Rho *
      (
          2.0 * mu * (epsilon(ux, uy, uz)' * epsilon(vx, vy, vz))
        + lameLambda * div(ux, uy, uz) * div(vx, vy, vz)
      )
)
+ on(2, ux = 0, uy = 0, uz = 0);

varf Mass([ux, uy, uz], [vx, vy, vz]) =
int3d(Th)(RhoMass * (ux * vx + uy * vy + uz * vz))
+ on(2, ux = 0, uy = 0, uz = 0);

int showeigsize = 1;
macro SolveFirstEigen(modeX, modeY, modeZ, lambda1)
{
    real eigenClockStart = clock();
    matrix OPm = EigenOP(Uh, Uh, tgv = 1e30);
    matrix Mm = Mass(Uh, Uh, tgv = 1.0);
    set(OPm, solver = sparsesolver);

    if (showeigsize == 1)
    {
        cout << "EigenOP matrix size = " << OPm.n << " x " << OPm.m << endl;
        object << "EigenOP matrix size = " << OPm.n << " x " << OPm.m << endl;
        showeigsize = 0;
    }

    real[int] eigVal(nev);
    Uh[int] [eigX, eigY, eigZ](nev);

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

    int eigCount = min(nconv, nev);

    lambda1 = eigVal[0];
    lambda2 = 0.0;
    lambda3 = 0.0;
    lambda4 = 0.0;
    lambda5 = 0.0;
    lambda6 = 0.0;
    gap12 = 1e30;
    gap12Ratio = 1e30;

    real denomGap = abs(lambda1);
    if (denomGap < 1e-30)
    {
        denomGap = 1e-30;
    }
    if (eigCount > 1)
    {
        lambda2 = eigVal[1];
        gap12 = lambda2 - lambda1;
        gap12Ratio = gap12 / denomGap;
    }
    if (eigCount > 2)
    {
        lambda3 = eigVal[2];
    }
    if (eigCount > 3)
    {
        lambda4 = eigVal[3];
    }
    if (eigCount > 4)
    {
        lambda5 = eigVal[4];
    }
    if (eigCount > 5)
    {
        lambda6 = eigVal[5];
    }

    cout << "Eigenvalues: lambda1 = " << lambda1
         << ", lambda2 = " << lambda2
         << ", lambda3 = " << lambda3
         << ", lambda4 = " << lambda4
         << ", lambda5 = " << lambda5
         << ", lambda6 = " << lambda6 << endl;
    cout << "converged eigenpairs used = " << eigCount << endl;
    cout << "gap12 = " << gap12
         << ", gap12/lambda1 = " << gap12Ratio << endl;

    object << "Eigenvalues: lambda1 = " << lambda1
           << ", lambda2 = " << lambda2
           << ", lambda3 = " << lambda3
           << ", lambda4 = " << lambda4
           << ", lambda5 = " << lambda5
           << ", lambda6 = " << lambda6 << endl;
    object << "converged eigenpairs used = " << eigCount << endl;
    object << "gap12 = " << gap12
           << ", gap12/lambda1 = " << gap12Ratio << endl;

    if ((nconv > 1) && (gap12Ratio < modeGapTol))
    {
        cout << "Warning: possible mode switching or nearly multiple first eigenvalue." << endl;
        object << "Warning: possible mode switching or nearly multiple first eigenvalue." << endl;
    }

    modeX[] = eigX[0][];

    real massNorm = sqrt(int3d(Th)(RhoMass * (modeX^2 + modeY^2 + modeZ^2)));
    if (massNorm > 1e-30)
    {
        modeX[] = modeX[] / massNorm;
    }

    lastEigenSolveTime = clock() - eigenClockStart;
    cout << "EigenOP Success!" << endl;
    cout << "Eigen solve time = " << lastEigenSolveTime << " s" << endl;
    object << "Eigen solve time = " << lastEigenSolveTime << " s" << endl;
} //
~~~

## 8. 第一特征值的形状导数

假设第一特征值 \(\lambda_1\) 是简单特征值，对应振型为 \(u_1\)，并满足归一化条件

\[
m_{\rho_M}(u_1,u_1)=1.
\]

特征值方程为

\[
a_{\rho_K}(u_1,v)=\lambda_1m_{\rho_M}(u_1,v).
\]

对密度求导时，由于广义特征值问题是自伴随的，振型导数项通过特征值方程抵消，得到

\[
d\lambda_1
=a_{\rho_K}'(u_1,u_1)-\lambda_1m_{\rho_M}'(u_1,u_1).
\]

其中

\[
a_{\rho_K}'(u_1,u_1)
=
\int_D
\left(
2\mu\varepsilon(u_1):\varepsilon(u_1)
+
\lambda_L(\operatorname{div}u_1)^2
\right)\delta\rho_K\,dx,
\]

\[
m_{\rho_M}'(u_1,u_1)
=
\int_D |u_1|^2\delta\rho_M\,dx.
\]

在可设计区域中

\[
\delta\rho_K=(1-\mathrm{weak})H_\gamma'(\phi)\delta\phi,
\qquad
\delta\rho_M=(1-\mathrm{weakM})H_\gamma'(\phi)\delta\phi.
\]

因此可以写成

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

水平集输运方程为

\[
\phi_t+\theta\cdot\nabla\phi=0.
\]

由于小质量块区域固定不变，在可设计区域有

\[
\rho_{K,t}
=-(1-\chi_{D_{\mathrm{block}}})H_\gamma'(\phi)\theta\cdot\nabla\phi
=-\theta\cdot\nabla\rho_{K,\mathrm{design}}.
\]

因此

\[
d\lambda_1(\theta)
=-\int_D S(u_1,\lambda_1)\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

令 \(J(\phi)=-\lambda_1(\phi)\)，则

\[
dJ(\theta)
=\int_D S(u_1,\lambda_1)\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

## 9. 保体积速度场

体积取为刚度密度对应的实体体积：

\[
\operatorname{Vol}(\phi)=\int_D\rho_K(\phi)\,dx.
\]

沿水平集输运方程演化时，小质量块区域不参与优化，所以

\[
\frac{d}{dt}\operatorname{Vol}(\phi)
=-\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

只要速度场满足

\[
\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0,
\]

就得到连续层面的体积守恒。

求 \(\theta\) 和约束变量 \(p\)，使得

\[
\int_D \nabla\theta:\nabla w\,dx
+
\int_D \theta\cdot w\,dx
+p\int_D w\cdot\nabla\rho_{K,\mathrm{design}}\,dx
=-\int_D S(u_1,\lambda_1)w\cdot\nabla\rho_{K,\mathrm{design}}\,dx,
\]

并满足

\[
\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0.
\]

为避免显式组装鞍点矩阵，程序采用 Schur 补形式。记

\[
A_\theta z=b_\rho,
\qquad
A_\theta y=c_\rho.
\]

则

\[
\theta=z-py,
\qquad
p=\frac{c_\rho^Tz}{c_\rho^Ty}.
\]

速度场方程是对称正定的向量型 Helmholtz 方程。代码对 `Atheta` 设置 `solver = sparsesolver`；由于前面加载了 `MUMPS_seq`，这里实际使用串行 MUMPS 后端做稀疏直接求解。对于当前由 `hMeshTarget = 0.09` 控制的三维网格，该设置通常比整体组装大鞍点系统更稳。

~~~freefem
Thetah [thetax, thetay, thetaz], [Wx, Wy, Wz];

varf VelocityA([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
int3d(Th)(
        Grad(thetax, thetay, thetaz)' * Grad(Wx, Wy, Wz)
      + [thetax, thetay, thetaz]' * [Wx, Wy, Wz]
)
+ on(2, thetax = 0)
+ on(3, thetax = 0, thetay = 0, thetaz = 0)
+ on(11, thetax = 0)
+ on(12, thetay = 0)
+ on(13, thetaz = 0);

matrix Atheta = VelocityA(Thetah, Thetah);
set(Atheta, solver = sparsesolver);

cout << "Velocity matrix size = " << Atheta.n << " x " << Atheta.m << endl;
object << "Velocity matrix size = " << Atheta.n << " x " << Atheta.m << endl;

macro VelocityPDE(thetax, thetay, thetaz, modex, modey, modez, lambdaone, phiLS)
{
    real velocityClockStart = clock();
    varf Velocityc([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
    int3d(Th)(
        [Wx, Wy, Wz]' * gradRhoDesign(phiLS)
    )
    + on(2, thetax = 0)
    + on(3, thetax = 0, thetay = 0, thetaz = 0)
    + on(11, thetax = 0)
    + on(12, thetay = 0)
    + on(13, thetaz = 0);

    varf Velocityb([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
    - int3d(Th)(
        (
            2.0 * mu * (epsilon(modex, modey, modez)' * epsilon(modex, modey, modez))
          + lameLambda * div(modex, modey, modez) * div(modex, modey, modez)
          - lambdaone * (1.0 - weakM) / (1.0 - weak) * (modex^2 + modey^2 + modez^2)
        )
        * ([Wx, Wy, Wz]' * gradRhoDesign(phiLS))
    )
    + on(2, thetax = 0)
    + on(3, thetax = 0, thetay = 0, thetaz = 0)
    + on(11, thetax = 0)
    + on(12, thetay = 0)
    + on(13, thetaz = 0);

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
    lastVelocitySolveTime = clock() - velocityClockStart;
    cout << "Velocity Success!" << endl;
    cout << "Velocity solve time = " << lastVelocitySolveTime << " s" << endl;
    object << "Velocity solve time = " << lastVelocitySolveTime << " s" << endl;
} //
~~~

## 10. 监测量

程序中主要记录以下量：

1. 第一特征值 \(\lambda_1\)；
2. 当前体积 \(\operatorname{Vol}(\phi)=\int_D\rho_K\,dx\)；
3. 速度场对应的体积变化率

\[
d\operatorname{Vol}(\theta)=\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx;
\]

4. 速度场对应的第一特征值变化率

\[
d\lambda_1(\theta)
=-\int_D S(u_1,\lambda_1)\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx.
\]

若速度场方程求解正确，则 `dVolTheta` 应接近零。若 `dLambdaTheta` 为正，则说明当前速度场方向会使第一特征值上升。

~~~freefem
macro GetVol(Rho)
(
    int3d(Th)(Rho)
) //

macro GetdVol(thetax, thetay, thetaz, phiLS)
(
    int3d(Th)([thetax, thetay, thetaz]' * gradRhoDesign(phiLS))
) //

macro GetdLambda(modeX, modeY, modeZ, thetax, thetay, thetaz, lambda1, phiLS)
(
    - int3d(Th)(
        (
            2.0 * mu * (epsilon(modeX, modeY, modeZ)' * epsilon(modeX, modeY, modeZ))
          + lameLambda * div(modeX, modeY, modeZ) * div(modeX, modeY, modeZ)
          - lambda1 * (1.0 - weakM) / (1.0 - weak) * (modeX^2 + modeY^2 + modeZ^2)
        )
        * ([thetax, thetay, thetaz]' * gradRhoDesign(phiLS))
      )
) //
~~~

## 11. 三维 dual-P^1 DG 水平集输运方程

水平集函数通过三维输运方程更新：

\[
\phi_t+\theta\cdot\nabla\phi=0.
\]

这里 \(\phi\) 使用间断有限元空间 `Phih=P1dc`。更新时采用隐式迎风 dual-P^1 DG 格式。单步离散为

\[
\int_D \frac{\phi^{n+1}}{\Delta t}v\,dx
+
\int_D (\theta\cdot\nabla\phi^{n+1})v\,dx
-
\sum_{F}\int_F
v\left(\frac{|\theta\cdot n|}{2}-\frac{\theta\cdot n}{2}\right)[\phi^{n+1}]\,ds
=
\int_D\frac{\phi^n}{\Delta t}v\,dx.
\]

代码中按照你给出的三维 dual-P^1 DG 写法使用 `int3d(Th)` 和 `intallfaces(Th)`，并把矩阵转换为 PETSc `Mat` 后求解。这里每次形状更新把 `Tsum` 分成 5 个 dual-P^1 DG 小步：

\[
\Delta t=\frac{Tsum}{5}.
\]

~~~freefem
macro NormalTheta(thetax, thetay, thetaz)
(
    N.x * thetax + N.y * thetay + N.z * thetaz
) //

macro LFDG(phiOutput, phiInput, thetax, thetay, thetaz, Tsum)
{
    real t = 0.0;
    real dt = Tsum / 5.0;

    Phih u, v;

    varf Adual(u, v) =
    int3d(Th)((u / dt + ([thetax, thetay, thetaz]' * grad(u))) * v)
    - intallfaces(Th)(
          (nElementonB - 1.0) * v
        * (0.5 * abs(NormalTheta(thetax, thetay, thetaz))
           - NormalTheta(thetax, thetay, thetaz) / 2.0)
        * jump(u)
      );

    varf DGRHS(u, v) = int3d(Th)((u / dt) * v);

    matrix Bm = DGRHS(Phih, Phih);
    matrix Am = Adual(Phih, Phih);
    Mat B(Bm), A(Am);

    for (t = 0.0; t < Tsum; t += dt)
    {
        v[] = B * phiInput[];
        phiOutput[] = A^-1 * v[];
        phiInput = phiOutput;
    }

    cout << "DG Success!" << endl;
} //
~~~

## 12. 重新初始化

重新初始化使用 FreeFEM 的 `distance`。输入可以是当前的 dual-P^1 DG 水平集函数，输出数组必须属于连续的 `P1h` 空间。因此代码中先定义

\[
\texttt{P1h phiRe},
\]

再执行

\[
\texttt{distance(Th, phiInput, phiRe[])}.
\]

最后把重新初始化结果赋值回 dual-P^1 DG 空间。

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
- 当前体积 \(\int_D\rho_K\,dx\)；
- 相比初始体积的体积变化；
- 当前速度场对应的体积变化率；
- 当前速度场对应的第一特征值变化率；
- 第二特征值以及第一、第二特征值相对间隔。

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
            << ";   diff from initial = " << vol[i] - vol[0]
            << ";   step diff = " << vol[i] - vol[i-1] << endl;
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

1. 已知当前 \(\phi^n\)，计算刚度密度和质量密度：

\[
\rho_K^n
=\chi_{D_{\mathrm{block}}}
+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weak}}(\phi^n),
\]

\[
\rho_M^n
=\mathrm{blockMass}\,\chi_{D_{\mathrm{block}}}
+(1-\chi_{D_{\mathrm{block}}})H_{\gamma,\mathrm{weakM}}(\phi^n).
\]

2. 解弹性结构第一特征值问题：

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

4. 用三维 dual-P^1 DG 格式求解水平集输运方程：

\[
\phi_t+\theta^n\cdot\nabla\phi=0.
\]

5. 每 20 次迭代使用 `distance` 重新初始化，使 \(\phi\) 保持接近有向距离函数。

~~~freefem
Phih phitest;
real[int] EigenValue1(iterMax + 1);
real[int] EigenValue2(iterMax + 1);
real[int] EigenValue3(iterMax + 1);
real[int] EigenValue4(iterMax + 1);
real[int] EigenValue5(iterMax + 1);
real[int] EigenValue6(iterMax + 1);
real[int] Gap12Ratio(iterMax + 1);
real[int] Vol(iterMax + 1);
real dVolTheta, dLambdaTheta;

int iter = 0;

SolveFirstEigen(modeX, modeY, modeZ, lambda1);
VelocityPDE(thetax, thetay, thetaz, modeX, modeY, modeZ, lambda1, phi);

EigenValue1[iter] = lambda1;
EigenValue2[iter] = lambda2;
EigenValue3[iter] = lambda3;
EigenValue4[iter] = lambda4;
EigenValue5[iter] = lambda5;
EigenValue6[iter] = lambda6;
Gap12Ratio[iter] = gap12Ratio;
Vol[iter] = GetVol(Rho);
dVolTheta = GetdVol(thetax, thetay, thetaz, phi);
dLambdaTheta = GetdLambda(modeX, modeY, modeZ,
                           thetax, thetay, thetaz,
                           lambda1, phi);
OutPutData(iter, EigenValue1, Vol, dVolTheta, dLambdaTheta, object);
plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter, wait = false);


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
    LFDG(phi, phitest, thetax, thetay, thetaz, Titer);

    if ((iter + 1) % 20 == 0)
    {
        RedistanceConvert(phi);
        cout << "Redistance End!" << endl;
    }

    Rho = RhoWithBlock(phi);
    RhoMass = RhoMassWithBlock(phi);

    SolveFirstEigen(modeX, modeY, modeZ, lambda1);
    VelocityPDE(thetax, thetay, thetaz, modeX, modeY, modeZ, lambda1, phi);

    iter = iter + 1;

    EigenValue1[iter] = lambda1;
    EigenValue2[iter] = lambda2;
    EigenValue3[iter] = lambda3;
    EigenValue4[iter] = lambda4;
    EigenValue5[iter] = lambda5;
    EigenValue6[iter] = lambda6;
    Gap12Ratio[iter] = gap12Ratio;
    Vol[iter] = GetVol(Rho);
    dVolTheta = GetdVol(thetax, thetay, thetaz, phi);
    dLambdaTheta = GetdLambda(modeX, modeY, modeZ,
                               thetax, thetay, thetaz,
                               lambda1, phi);

    cout << "Titer = " << Titer << endl;
    object << "Titer = " << Titer << endl;
    OutPutData(iter, EigenValue1, Vol, dVolTheta, dLambdaTheta, object);

    
    plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter, wait = false);
    
    if (iter % 10 == 0)
    {
        string filename = "./figure/Rho3D_" + iter + ".vtu";
        savevtk(filename, Th, Rho);

        string filenamePhi = "./figure/Phi3D_" + iter + ".vtu";
        savevtk(filenamePhi, Th, phi);
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

最终输出 Matlab 数据，包含每一步的体积、前六个特征值以及第一、第二特征值间隔比。

其中

\[
\texttt{Volume(i)}
\]

对应第 \(i\) 步的体积，

\[
\texttt{lambda1(i)}
\]

对应第 \(i\) 步的第一特征值，`gap12ratio(i)` 用于检测模态切换或近重特征值。

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
    OutMatlab << "lambda5(" << i+1 << ") = " << EigenValue5[i] << "; ";
    OutMatlab << "lambda6(" << i+1 << ") = " << EigenValue6[i] << "; ";
    OutMatlab << "gap12ratio(" << i+1 << ") = " << Gap12Ratio[i] << ";" << endl;
}
~~~

## 16. 模型总结

这个程序实际求解的是三维第一特征值最大化问题：

\[
\max_\phi \lambda_1(\phi),
\]

其中 \(\lambda_1(\phi)\) 由三维广义特征值问题给出：

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
-
\int_D S(u_1,\lambda_1)w\cdot\nabla\rho_{K,\mathrm{design}}\,dx,
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

因为速度场满足

\[
\int_D\theta\cdot\nabla\rho_{K,\mathrm{design}}\,dx=0,
\]

所以体积在连续层面满足

\[
\frac{d}{dt}\int_D\rho_K\,dx=0.
\]

这就是该程序不需要额外添加任意体积乘子的原因。体积约束已经包含在速度场方程中。

## 17. dual-P^1 DG 求解耗时、求解器对比与 T 策略

当前程序额外输出 dual-P^1 DG 水平集输运方程的求解耗时，并保留 Eigen、Velocity 和单步总耗时，便于判断每一步的主要耗时来源。正式默认参数已经改为第三种 T 策略：

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `iterMax` | `300` | 正式优化步数 |
| `T0` | `0.2` | T 的上界 |
| `TStrategy` | `3` | 使用 \(T=\min(T_0,\ell h/\|\theta\|_\infty)\) |
| `TScale` | `0.2` | 限制单步最大界面位移约为 `0.2h` |
| `TDecayPeriod` | `5` | 仅策略 2 使用 |
| `dual-P^1 DGSolverMode` | `1` | dual-P^1 DG 方程使用 GMRES；`0` 为 MUMPS_seq/direct |
| `RedistancePeriod` | `20` | 每隔 20 步重初始化；重初始化步允许出现体积跳变，重点检查非重初始化步 |
| `InitialVolumeRatio` | `0.5` | 初始形状目标体积占整体区域体积的比例 |
| `InitialVolumeTol` | `1e-5` | 初始体积校准容差预留参数 |

默认正式运行命令为：

~~~powershell
FreeFem++ -v 0 Main.edp
~~~

也可以显式写出同一组参数：

~~~powershell
FreeFem++ -v 0 Main.edp -iterMax 300 -T0 0.2 -TStrategy 3 -TScale 0.2 -dual-P^1 DGSolverMode 1 -RedistancePeriod 20 -InitialVolumeRatio 0.5
~~~

其中 `dual-P^1 DGSolverMode = 0` 表示 dual-P^1 DG 方程使用 `MUMPS_seq/sparsesolver` 直接法，`dual-P^1 DGSolverMode = 1` 表示 dual-P^1 DG 方程使用 `GMRES`。在 `hMeshTarget = 0.09`、单步测试中，dual-P^1 DG 直接法约 41 秒，GMRES 约 16 秒；整体单步从约 135 秒降到约 110 秒。Eigen 和 Velocity 仍保持当前 `sparsesolver`，所以该对比主要反映 dual-P^1 DG 方程的直接法与迭代法差异。

程序还输出以下诊断量：

- `theta L2 norm` 和 `theta max norm`：判断速度场尺度；
- `T*max|theta|/h`：判断单步界面位移相对网格尺寸的大小；
- `Phi L2 change`：判断 dual-P^1 DG 更新后的水平集实际变化；
- `Vol step diff`：判断每一步体积变化是否稳定在目标范围内。

两步短测试结果显示，`T0 = 0.2`、恒定 T、dual-P^1 DG 使用 GMRES 时：

- 第 1 步 `T*max|theta|/h = 0.1604`，`Vol step diff = -1.7531e-4`；
- 第 2 步 `T*max|theta|/h = 0.142865`，`Vol step diff = -1.12673e-4`；
- 两步体积单步变化均明显小于 `0.002`，`lambda1` 连续上升。

因此推荐的保守正式策略是

~~~powershell
FreeFem++ -v 0 Main.edp -iterMax 300 -T0 0.2 -TStrategy 3 -TScale 0.2 -dual-P^1 DGSolverMode 1 -RedistancePeriod 20 -InitialVolumeRatio 0.5
~~~

这里 `TStrategy = 3` 使用

\[
T=\min(T_0,\ell h/\|\theta\|_\infty),
\]

即使后续速度场突然变大，也把单步最大界面位移限制在约 `0.2h` 以内。`T0 = 0.6` 的单步体积变化仍小于 `0.002`，但第一步已经达到 `T*max|theta|/h = 0.481199`，不建议作为默认连续迭代参数，避免第二步以后过大更新导致不稳定。

`distance` 重初始化可能改变水平集零等值面附近的平滑密度积分，导致重初始化步出现明显体积跳变。因此体积稳定性主要看非重初始化步；若重初始化后速度场或体积变化过大，可通过第三种 T 策略自动限制步长，必要时再调小 `TScale`，但不应突然大幅降低 `T0`。

如果形状变化仍不明显，应优先检查速度场方程与线弹性特征值敏感度，而不是先怀疑 dual-P^1 DG PDE：

- 当前 `dVol(theta)` 在测试中为 `1e-16` 量级，说明速度场的体积投影约束基本生效；
- `dLambda(theta)` 为正且 `lambda1` 连续上升，说明速度方向与目标函数一阶变化方向一致；
- 若 `T*max|theta|/h` 已达到 `0.15` 到 `0.2` 但几何变化仍弱，应检查速度场 RHS 的尺度、特征模态归一化、质量密度系数 \((1-\mathrm{weakM})/(1-\mathrm{weak})\)，以及保体积投影中 `cRho`、`bRho` 的符号约定。
