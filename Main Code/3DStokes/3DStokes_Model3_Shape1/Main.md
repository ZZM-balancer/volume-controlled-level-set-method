# 带详细说明的 FreeFEM 代码：Brinkman-Stokes 保体积水平集拓扑优化

这个文件的组织方式是：**普通 Markdown 说明文字 + FreeFEM 代码块**。

- Markdown 说明部分用于解释数学模型、变量含义、线性系统、求解器设置和每一步计算的意义。
- 每个 `~~~freefem ... ~~~` 代码块内部尽量保持为纯 FreeFEM 代码。
- 若需要运行程序，可以把所有 FreeFEM 代码块按本文顺序连起来，得到完整 `.edp` 文件。
- 本版本保持你的要求：固定迭代时间 `Titer`，不做下降性检测，不做步长回退。
- 本版本将原先的 `grad(Rho)` 替换为更接近连续公式的

\[ \nabla \rho = \nabla H_\gamma(\phi) \approx H_\gamma'(\phi)\nabla \phi . \]

这样做的原因是：先把 `Rho = H_gamma(phi)` 插值到 `P1dc3d`，再取 `grad(Rho)`，实际上是在对插值后的密度函数求梯度；而形状导数中的连续表达更自然的是

\[ H_\gamma'(\phi)\nabla \phi . \]

因此，下降速度场、体积约束、方向导数输出中都统一使用 `gradRhoPhi(phi)`。

---

## 0. 模型总览

设计区域为

\[ D\subset \mathbb{R}^3. \]

水平集函数为

\[ \phi:D\to \mathbb{R}, \]

形状区域取为

\[ \Omega=\{x\in D:\phi(x)<0\}. \]

密度变量为

\[ \rho=H_\gamma(\phi), \]

其中

\[ H_\gamma(\phi)= \begin{cases} 1, & \phi< -\gamma,\\[4pt] \dfrac12-\dfrac{\phi}{2\gamma} -\dfrac{1}{2\pi}\sin\left(\dfrac{\pi\phi}{\gamma}\right), & -\gamma\leq \phi\leq \gamma,\\[8pt] 0, & \phi>\gamma. \end{cases} \]

所以

\[ \rho\approx 1 \quad\Longleftrightarrow\quad \phi<0, \]

\[ \rho\approx 0 \quad\Longleftrightarrow\quad \phi>0. \]

Brinkman 惩罚系数为

\[ \alpha(\rho) = \overline{\alpha} + (\underline{\alpha}-\overline{\alpha}) \frac{(1+q)\rho}{q+\rho}, \]

其导数为

\[ \alpha'(\rho) = (\underline{\alpha}-\overline{\alpha}) \frac{q(1+q)}{(q+\rho)^2}. \]

状态方程为 Brinkman-Stokes 方程：给定 \(\rho\)，求 \((u,p)\)，使得

\[ \begin{cases} -\Delta u+\alpha(\rho)u+\nabla p=f, & \text{in }D,\\[3pt] \operatorname{div}u=0, & \text{in }D,\\[3pt] u=g, & \text{on }\Gamma_{\mathrm{in}},\\[3pt] u=0, & \text{on }\Gamma_{\mathrm{wall}}. \end{cases} \]

目标函数为

\[ J(\rho,u) = \frac12\int_D \alpha(\rho)|u|^2\,dx + \frac12\int_D |Du|^2\,dx - \int_D f\cdot u\,dx. \]

在当前代码中 \(f=0\)，因此数值上就是

\[ J(\rho,u) = \frac12\int_D \alpha(\rho)|u|^2\,dx + \frac12\int_D |Du|^2\,dx. \]

体积为

\[ \operatorname{Vol}(\phi)=\int_D H_\gamma(\phi)\,dx. \]

下降速度场 \(\theta\) 通过带体积约束的辅助问题求解。记

\[ a(\theta,w) = \int_D D\theta:Dw\,dx + \int_D \theta\cdot w\,dx, \]

\[ c(w)=\int_D w\cdot H_\gamma'(\phi)\nabla\phi\,dx, \]

\[ \ell(w) = \int_D \frac12\alpha'(\rho)|u|^2 w\cdot H_\gamma'(\phi)\nabla\phi\,dx. \]

求

\[ (\theta,\lambda)\in [H_0^1(D)]^3\times\mathbb{R}, \]

使得

\[ a(\theta,w)+\lambda c(w)=\ell(w), \qquad \forall w\in [H_0^1(D)]^3, \]

并且

\[ c(\theta)=0. \]

这个约束即

\[ \int_D \theta\cdot H_\gamma'(\phi)\nabla\phi\,dx=0, \]

用于近似体积守恒方向。

水平集更新方程为

\[ \frac{\partial \phi}{\partial t} + \theta\cdot\nabla\phi =0. \]

每次优化迭代中，时间长度固定为

\[ T_{\mathrm{iter}}. \]

本代码不做任何下降性检测或步长回退，因此每次都直接接受由固定 `Titer` 得到的新 \(\phi\)。

---

## 1. 加载模块与全局参数

这一段负责加载 FreeFEM 模块并设置全局参数。

### 主要参数含义

| 变量 | 含义 |
|---|---|
| `Titer` | 每次优化迭代中水平集演化的总时间 |
| `iterMax` | 最大优化迭代次数 |
| `levelSubSteps` | 每个 `Titer` 内部再分成几个 DG 子步；为了速度默认设为 1 |
| `redistEvery` | 每隔多少次优化迭代做一次重初始化 |
| `saveEvery` | 每隔多少次迭代保存一次 `.vtu` 文件 |
| `alphaLow` | \(\underline{\alpha}\)，低阻力区域惩罚系数 |
| `alphaUp` | \(\overline{\alpha}\)，高阻力区域惩罚系数 |
| `qq` | Brinkman 插值函数中的参数 \(q\) |
| `gamma` | 平滑 Heaviside 函数的半宽度 |
| `useMonitor` | 是否打开 PETSc 迭代残差输出；`0` 表示关闭，`1` 表示打开 |

`levelSubSteps = 1` 是为了尽量把每次迭代时间压低。若后续发现界面演化不够平滑，可以改为 `2`，但计算时间会增加。

~~~freefem
load "TetGen"
load "gmsh"
load "distance"
load "iovtk"
load "PETSc"

verbosity = 0;

real Titer = 0.001;
int iterMax = 300;

int levelSubSteps = 5;
int redistEvery = 20;
int saveEvery = 10;

real alphaLow = 2.5e-4;
real alphaUp  = 2.5e+4;
real qq = 0.1;
real gamma = 0.005;

ofstream object("./Output/Main.obj");

string ModelShapeName  = "figure/InitialShape.vtu";
string ModelResultName = "figure/FinalShape.vtu";
string ModelMatLabName = "matlab/Model3_Shape1.m";

int useMonitor = 0;
~~~

---

## 2. 基本宏定义

这一段定义密度函数、惩罚系数、导数和常用微分算子。

### 2.1 平滑 Heaviside 函数

`HeaviSide(phi)` 对应

\[ H_\gamma(\phi)= \begin{cases} 1, & \phi< -\gamma,\\[3pt] \dfrac12-\dfrac{\phi}{2\gamma} -\dfrac{1}{2\pi}\sin\left(\dfrac{\pi\phi}{\gamma}\right), & |\phi|\leq \gamma,\\[8pt] 0, & \phi>\gamma. \end{cases} \]

### 2.2 `dHeaviSide(phi)`

这是本版代码相对原代码的关键修改。由上式求导得

\[ H_\gamma'(\phi) = -\frac{1}{2\gamma} \left(1+\cos\frac{\pi\phi}{\gamma}\right), \qquad |\phi|\leq \gamma, \]

并且

\[ H_\gamma'(\phi)=0, \qquad |\phi|>\gamma. \]

因此

\[ \nabla \rho = \nabla H_\gamma(\phi) = H_\gamma'(\phi)\nabla\phi. \]

代码中使用 `gradRhoPhi(phi)` 表示这个量，即

\[ \texttt{gradRhoPhi(phi)} = H_\gamma'(\phi) \begin{bmatrix} \partial_x\phi\\ \partial_y\phi\\ \partial_z\phi \end{bmatrix}. \]

### 2.3 Brinkman 系数

`alpha(Rho)` 对应

\[ \alpha(\rho)= \overline{\alpha} + (\underline{\alpha}-\overline{\alpha}) \frac{(1+q)\rho}{q+\rho}. \]

`dalpha(Rho)` 对应

\[ \alpha'(\rho)= (\underline{\alpha}-\overline{\alpha}) \frac{q(1+q)}{(q+\rho)^2}. \]

由于 \(\underline{\alpha}<\overline{\alpha}\)，所以通常

\[ \alpha'(\rho)<0. \]

### 2.4 微分算子

- `grad(u)` 表示标量函数的梯度。
- `Grad(ux,uy,uz)` 表示三维向量场每个分量的梯度拼接。
- `div(ux,uy,uz)` 表示散度。

~~~freefem
macro HeaviSide(phi)
(
    (phi > gamma) * 0.0
  + (phi < -gamma) * 1.0
  + ((phi >= -gamma) && (phi <= gamma))
    * (0.5 - 0.5*phi/gamma - (0.5/pi)*sin(pi*phi/gamma))
)//

macro dHeaviSide(phi)
(
    ((phi >= -gamma) && (phi <= gamma))
    * ( -0.5/gamma * (1.0 + cos(pi*phi/gamma)) )
)//

macro gradRhoPhi(phi)
[
    dHeaviSide(phi) * dx(phi),
    dHeaviSide(phi) * dy(phi),
    dHeaviSide(phi) * dz(phi)
]//

macro alpha(Rho)
(
    alphaUp + (alphaLow - alphaUp) * (qq + 1.0) * Rho / (qq + Rho)
)//

macro dalpha(Rho)
(
    (alphaLow - alphaUp) * (qq + 1.0) * qq / ((qq + Rho)*(qq + Rho))
)//

macro grad(u)
[
    dx(u), dy(u), dz(u)
]//

macro Grad(ux, uy, uz)
[
    dx(ux), dy(ux), dz(ux),
    dx(uy), dy(uy), dz(uy),
    dx(uz), dy(uz), dz(uz)
]//

macro div(ux, uy, uz)
(
    dx(ux) + dy(uy) + dz(uz)
)//
~~~

---

## 3. 网格、边界标签与有限元空间

这一段读入三维网格并重新设置边界标签。

边界标签由 `Labelfunc` 给出：

\[ 1: y=0, \qquad 2: y=2, \qquad 3: \text{其余边界}. \]

因此可以理解为

\[ \Gamma_{\mathrm{in}}=\{y=0\}, \]

\[ \Gamma_{\mathrm{out}}=\{y=2\}, \]

\[ \Gamma_{\mathrm{wall}}=\partial D\setminus (\Gamma_{\mathrm{in}}\cup \Gamma_{\mathrm{out}}). \]

### 有限元空间

| 空间 | FreeFEM 写法 | 数学含义 | 用途 |
|---|---|---|---|
| `Uh` | `[P2,P2,P2]` | \([P_2]^3\) | Stokes 速度 \(u\) |
| `P1h` | `P1` | \(P_1\) | 压力 \(p\)、距离函数中间变量 |
| `Sth` | `Uh * P1h` | \([P_2]^3\times P_1\) | Stokes 混合空间 |
| `Thetah` | `[P1,P1,P1]` | \([P_1]^3\) | 下降速度场 \(\theta\) |
| `Phih` | `P1dc3d` | 不连续 \(P_1\) | 水平集函数 \(\phi\)、密度 \(\rho\) |

根据你给出的规模，当前问题大约为：

\[ \dim(S_h)=209810, \qquad \dim(\Theta_h)=27150. \]

因为

\[ \Theta_h=[P_1]^3, \]

所以 P1 节点数约为

\[ N_{P1}=27150/3=9050. \]

又因为

\[ S_h=[P_2]^3\times P_1, \]

所以

\[ 209810=3N_{P2}+9050, \]

即

\[ N_{P2}=66920. \]

`Phih=P1dc3d` 每个四面体单元有 4 个局部自由度，不同单元不共享自由度，因此

\[ \dim(\Phi_h)=4N_T. \]

对于你这个网格，预计

\[ \dim(\Phi_h)\approx 1.9\times 10^5. \]

也就是 level-set DG 线性系统规模大约在 `190000--200000`。

~~~freefem
mesh3 Th = gmshload3("Th006_box.msh");

macro Labelfunc()
(
      3
    - 2 * (abs(y) <= 1e-6)
    - 1 * (abs(y - 2.0) <= 1e-6)
)//

Th = change(Th, flabel = Labelfunc);

fespace Uh(Th, [P2, P2, P2]);
fespace P1h(Th, P1);
fespace Sth = Uh * P1h;

fespace Thetah(Th, [P1, P1, P1]);
fespace Phih(Th, P1dc3d);

Uh [uTmpx, uTmpy, uTmpz];
P1h pTmp;
Thetah [thetaTmpx, thetaTmpy, thetaTmpz];
Phih phiTmp;

object << "Sth total dof = " << (uTmpx[].n + pTmp[].n) << endl;
object << "Velocity field dof = " << thetaTmpx[].n << endl;
object << "Level-set DG dof = " << phiTmp[].n << endl;
object.flush;

cout << "Sth total dof = " << (uTmpx[].n + pTmp[].n) << endl;
cout << "Velocity field dof = " << thetaTmpx[].n << endl;
cout << "Level-set DG dof = " << phiTmp[].n << endl;

P1h h1 = hTriangle;
object << "hmax = " << h1[].max << endl;
object << "hmin = " << h1[].min << endl;
object.flush;
~~~

---

## 4. 初始水平集函数与初始密度

初始函数取为

\[ \phi_0(x,y,z) = 0.24- \sqrt{(y-0.6)^2+(x-0.5)^2+(z-0.5)^2}. \]

因此初始正值区域大致是以

\[ (0.5,0.6,0.5) \]

为中心、半径为 \(0.24\) 的球。

代码随后调用 FreeFEM 的 `distance` 命令，把这个初始函数转成一个近似有向距离函数。然后设置

\[ \rho^0=H_\gamma(\phi^0). \]

最后保存初始形状到 `InitialShape.vtu`。

~~~freefem
func phi0 = 0.2 - sqrt((y - 0.5)^2 + (x - 0.5)^2 + (z - 0.5)^2);

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho;
phi = phiInit;
Rho = HeaviSide(phi);

plot(Rho, fill = 1, value = 1, cmm = "Rho_0", wait = false);
savevtk(ModelShapeName, Th, phi, Rho);
~~~

---

## 5. Brinkman-Stokes 状态方程

给定当前密度

\[ \rho=H_\gamma(\phi), \]

求速度和压力

\[ (u,p):D\to \mathbb{R}^3\times\mathbb{R}, \]

满足

\[ \begin{cases} -\Delta u+\alpha(\rho)u+\nabla p=f, & \text{in }D,\\[3pt] \operatorname{div}u=0, & \text{in }D,\\[3pt] u=g, & \text{on }\Gamma_{\mathrm{in}},\\[3pt] u=0, & \text{on }\Gamma_{\mathrm{wall}}. \end{cases} \]

当前代码中

\[ f=(0,0,0), \]

并且入口速度为

\[ g=(0,32x(1-x)z(1-z),0). \]

出口 \(\Gamma_{\mathrm{out}}\) 上没有强制 Dirichlet 条件，因此对应自然边界条件。

### 弱形式

求

\[ (u_h,p_h)\in [P_2]^3\times P_1, \]

满足边界条件，并且对任意测试函数 \((v_h,q_h)\) 有

\[ \int_D Du_h:Dv_h\,dx - \int_D p_h\operatorname{div}v_h\,dx - \int_D q_h\operatorname{div}u_h\,dx -10^{-8}\int_D p_hq_h\,dx + \int_D \alpha(\rho)u_h\cdot v_h\,dx = \int_D f\cdot v_h\,dx. \]

这里的小项

\[ -10^{-8}\int_D p_hq_h\,dx \]

相当于压力稳定/压力零均值替代项，可以防止压力块奇异。

### 矩阵结构

线性系统为

\[ \begin{bmatrix} A+M_{\alpha(\rho)} & B^T\\ B & -\varepsilon M_p \end{bmatrix} \begin{bmatrix} U\\P \end{bmatrix} = \begin{bmatrix} F\\0 \end{bmatrix}, \qquad \varepsilon=10^{-8}. \]

其中：

- \(A\) 来自 \(\int_D Du:Dv\,dx\)。
- \(M_{\alpha(\rho)}\) 来自 \(\int_D \alpha(\rho)u\cdot v\,dx\)。
- \(B\) 来自散度约束。
- \(-\varepsilon M_p\) 是压力稳定项。

### 求解器设置

当前默认使用 PETSc：

\[ \text{GMRES} + \text{fieldsplit} + \text{Schur complement}. \]

速度块默认使用 HYPRE/BoomerAMG：

- `-fieldsplit_velocity_pc_type hypre`
- `-fieldsplit_velocity_pc_hypre_type boomeramg`

压力块使用 Jacobi：

- `-fieldsplit_pressure_pc_type jacobi`

如果你的 PETSc 没有 HYPRE，可以把速度块改成 GAMG，例如把速度块预处理改为 `gamg`。

`SolveStokes` 每次迭代只重新组装依赖 \(\rho\) 的 Brinkman 质量项，然后与固定 Stokes 部分相加。这样比每次完整重新组装所有项更快。

~~~freefem
P1h p, q;
Uh [ux, uy, uz], [vx, vy, vz];

func gx = 0.0;
func gy = 48.0*x*(1.0 - x)*z*(1.0 - z);
func gz = 0.0;

func fx = 0.0;
func fy = 0.0;
func fz = 0.0;

varf StokesSteady(<[ux, uy, uz], [p]>, <[vx, vy, vz], [q]>) =
int3d(Th)(
      Grad(ux, uy, uz)' * Grad(vx, vy, vz)
    - p * div(vx, vy, vz)
    - q * div(ux, uy, uz)
    - 1e-8 * p * q
)
+ on(1, ux = gx, uy = gy, uz = gz)
+ on(3, ux = 0.0, uy = 0.0, uz = 0.0);

varf StokesRHS(<[ux, uy, uz], [p]>, <[vx, vy, vz], [q]>) =
int3d(Th)(
      fx * vx
    + fy * vy
    + fz * vz
)
+ on(1, ux = gx, uy = gy, uz = gz)
+ on(3, ux = 0.0, uy = 0.0, uz = 0.0);

matrix AStokes0 = StokesSteady(Sth, Sth);
real[int] bStokes = StokesRHS(0, Sth);

cout << "Stokes system size = " << AStokes0.n << endl;
object << "Stokes system size = " << AStokes0.n << endl;

Uh [tagx, tagy, tagz];
P1h tagp;
[tagx, tagy, tagz] = [1.0, 1.0, 1.0];
tagp = 2.0;

real[int] fieldTags = [tagx[], tagp[]];
string[int] fieldNames(2);
fieldNames[0] = "velocity";
fieldNames[1] = "pressure";

string spStokes =
      "-ksp_type gmres "
    + "-ksp_gmres_restart 80 "
    + "-ksp_rtol 1e-6 "
    + "-ksp_atol 1e-9 "
    + "-ksp_max_it 400 "
    + "-pc_type fieldsplit "
    + "-pc_fieldsplit_type schur "
    + "-pc_fieldsplit_schur_precondition selfp "
    + "-fieldsplit_velocity_ksp_type preonly "
    /*+ "-fieldsplit_velocity_pc_type hypre "
    + "-fieldsplit_velocity_pc_hypre_type boomeramg "*/
    + "-fieldsplit_pressure_ksp_type preonly "
    + "-fieldsplit_pressure_pc_type jacobi ";

if (useMonitor == 1)
{
    spStokes = spStokes + "-ksp_monitor_true_residual -ksp_converged_reason ";
}

macro SolveStokes(ux, uy, uz, Rho)
{
    varf StokesRho(<[ux, uy, uz], [p]>, <[vx, vy, vz], [q]>) =
    int3d(Th)(
        alpha(Rho) * ([ux, uy, uz]' * [vx, vy, vz])
    );

    matrix ARho = StokesRho(Sth, Sth);
    matrix AStokes = AStokes0 + ARho;

    Mat APetsc(AStokes);
    set(APetsc, sparams = spStokes, fields = fieldTags, names = fieldNames);

    real[int] sol = APetsc^-1 * bStokes;

    int nU = ux[].n;
    ux[] = sol(0:nU - 1);
    p[]  = sol(nU:sol.n - 1);
    cout << "Stokes Success! "<< endl;
}//
~~~

---

## 6. 下降速度场方程

这一段求解形状更新方向 \(\theta\)。

### 原始带约束问题

定义

\[ a(\theta,w) = \int_D D\theta:Dw\,dx + \int_D \theta\cdot w\,dx. \]

体积约束泛函为

\[ c(w) = \int_D w\cdot H_\gamma'(\phi)\nabla\phi\,dx. \]

右端项为

\[ \ell(w) = \int_D \frac12\alpha'(\rho)|u|^2 w\cdot H_\gamma'(\phi)\nabla\phi\,dx. \]

求

\[ (\theta,\lambda)\in [H_0^1(D)]^3\times\mathbb{R}, \]

使得

\[ a(\theta,w)+\lambda c(w)=\ell(w), \qquad \forall w\in [H_0^1(D)]^3, \]

并满足

\[ c(\theta)=0. \]

矩阵形式为

\[ \begin{bmatrix} A_\theta & c_\rho\\ c_\rho^T & 0 \end{bmatrix} \begin{bmatrix} \Theta\\ \lambda \end{bmatrix} = \begin{bmatrix} b_\rho\\0 \end{bmatrix}. \]

### 本代码采用的快速求法

因为 \(\lambda\) 只是一个标量，不必直接求解整个鞍点系统。代码中改为两次 SPD 求解：

\[ A_\theta y=b_\rho, \]

\[ A_\theta z=c_\rho. \]

然后由约束

\[ c_\rho^T\theta=0 \]

得到

\[ \lambda=\frac{c_\rho^T y}{c_\rho^T z}, \]

再令

\[ \theta=y-\lambda z. \]

这样做的好处是：

1. 不需要求解不定的鞍点系统；
2. `Atheta` 是固定矩阵，可以重复使用；
3. 对你给出的速度场规模 `27150`，两次 SPD 求解通常比一个鞍点系统更稳定。

### 与原代码的关键区别

原代码中大量使用

\[ \nabla\rho \]

并通过 `grad(Rho)` 计算。这里统一替换为

\[ H_\gamma'(\phi)\nabla\phi. \]

对应代码就是 `gradRhoPhi(phi)`。

~~~freefem
Thetah [thetax, thetay, thetaz], [Wx, Wy, Wz];

varf VelocityA([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
int3d(Th)(
      Grad(thetax, thetay, thetaz)' * Grad(Wx, Wy, Wz)
    + [thetax, thetay, thetaz]' * [Wx, Wy, Wz]
)
+ on(1, 2, 3, thetax = 0.0, thetay = 0.0, thetaz = 0.0);

matrix Atheta = VelocityA(Thetah, Thetah);
cout << "Velocity system size = " << Atheta.n << endl;
object << "Velocity system size = " << Atheta.n << endl;

Mat AthetaPetsc(Atheta);

string spTheta =
      "-ksp_type cg "
    + "-ksp_rtol 1e-8 "
    + "-ksp_atol 1e-10 "
    + "-ksp_max_it 300 "
    + "-pc_type hypre "
    + "-pc_hypre_type boomeramg ";

if (useMonitor == 1)
{
    spTheta = spTheta + "-ksp_monitor_true_residual -ksp_converged_reason ";
}

set(AthetaPetsc, sparams = spTheta);

macro SolveVelocity(thetax, thetay, thetaz, ux, uy, uz, Rho, phi)
{
    varf VelocityC([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
    int3d(Th)(
        [Wx, Wy, Wz]' * gradRhoPhi(phi)
    )
    + on(1, 2, 3, thetax = 0.0, thetay = 0.0, thetaz = 0.0);

    varf VelocityB([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
    int3d(Th)(
        dalpha(Rho)
        * 0.5 * ([ux, uy, uz]' * [ux, uy, uz])
        * ([Wx, Wy, Wz]' * gradRhoPhi(phi))
    )
    + on(1, 2, 3, thetax = 0.0, thetay = 0.0, thetaz = 0.0);

    real[int] cRho = VelocityC(0, Thetah);
    real[int] bRho = VelocityB(0, Thetah);

    real[int] ysol(bRho.n);
    real[int] zsol(cRho.n);

    ysol = AthetaPetsc^-1 * bRho;
    zsol = AthetaPetsc^-1 * cRho;

    real numerator = cRho' * ysol;
    real denominator = cRho' * zsol;
    real lambda = 0.0;

    if (abs(denominator) > 1e-30)
    {
        lambda = numerator / denominator;
    }

    real[int] thetaSol(ysol.n);
    for (int i = 0; i < ysol.n; i++)
    {
        thetaSol[i] = ysol[i] - lambda * zsol[i];
    }

    thetax[] = thetaSol;

    cout << "Velocity Success! "<< endl;
    object << "lambda = " << lambda << endl;
}//
~~~

---

## 7. 目标函数、体积和方向导数输出

这一段定义若干输出量，用于监测优化过程。

### 7.1 目标函数 `GetDrag`

代码计算

\[ J(\rho,u) = \frac12\int_D \alpha(\rho)|u|^2\,dx + \frac12\int_D |Du|^2\,dx - \int_D f\cdot u\,dx. \]

当前 \(f=0\)，所以第三项实际为零。

### 7.2 体积 `GetVol`

代码计算

\[ \operatorname{Vol}(\phi) = \int_D \rho\,dx = \int_D H_\gamma(\phi)\,dx. \]

这个值用于观察保体积效果。

### 7.3 方向导数 `GetdJ`

代码输出的是

\[ \operatorname{GetdJ}(\theta) = \int_D \left(-\frac12\alpha'(\rho)|u|^2\right) \theta\cdot H_\gamma'(\phi)\nabla\phi\,dx. \]

注意这里的符号与 `SolveVelocity` 中的右端项配套使用，主要用于检查当前方向的作用。

### 7.4 体积方向导数 `GetdV`

代码输出

\[ \operatorname{GetdV}(\theta) = \int_D \theta\cdot H_\gamma'(\phi)\nabla\phi\,dx. \]

因为下降速度场中加了约束，所以理想情况下应当有

\[ \operatorname{GetdV}(\theta)\approx 0. \]

如果这个量长期明显偏离 0，说明体积约束方向或者线性求解精度需要检查。

~~~freefem
macro GetDrag(ux, uy, uz, Rho)
(
    int3d(Th)(
        0.5 * alpha(Rho) * ([ux, uy, uz]' * [ux, uy, uz])
    )
  + int3d(Th)(
        0.5 * (Grad(ux, uy, uz)' * Grad(ux, uy, uz))
    )
  - int3d(Th)(
        fx * ux + fy * uy + fz * uz
    )
)//

macro GetVol(Rho)
(
    int3d(Th)(Rho)
)//

macro GetdJ(ux, uy, uz, thetax, thetay, thetaz, Rho, phi)
(
    int3d(Th)(
        (
            -dalpha(Rho)
            * 0.5 * ([ux, uy, uz]' * [ux, uy, uz])
        )
        * ([thetax, thetay, thetaz]' * gradRhoPhi(phi))
    )
)//

macro GetdV(thetax, thetay, thetaz, phi)
(
    int3d(Th)(
        [thetax, thetay, thetaz]' * gradRhoPhi(phi)
    )
)//
~~~

---

## 8. dual-P^1 DG 水平集演化方程

给定当前下降速度场 \(\theta\)，水平集函数通过对流方程更新：

\[ \frac{\partial\phi}{\partial t} + \theta\cdot\nabla\phi=0. \]

每次优化迭代演化固定时间

\[ T_{\mathrm{iter}}. \]

若 `levelSubSteps = m`，则每个 DG 子步时间为

\[ \Delta t=\frac{T_{\mathrm{iter}}}{m}. \]

当前为了控制每次迭代时间，默认

\[ m=1. \]

### dual-P^1 DG 离散形式

在一个子步中，求

\[ \phi_h^{k+1}\in P_1^{\mathrm{dc}}, \]

使得对任意测试函数 \(v_h\in P_1^{\mathrm{dc}}\)，有

\[ \int_D \frac{\phi_h^{k+1}}{\Delta t}v_h\,dx + \int_D (\theta\cdot\nabla\phi_h^{k+1})v_h\,dx - \sum_{e\in\mathcal{F}_h^i} \int_e v_h \frac{|\theta\cdot n|-\theta\cdot n}{2} [\phi_h^{k+1}]\,ds = \int_D \frac{\phi_h^k}{\Delta t}v_h\,dx. \]

其中：

- 第一项是隐式 Euler 时间项；
- 第二项是单元内部对流项；
- 第三项是内部面上的 upwind DG 数值通量项；
- 右端项来自上一子步的 \(\phi_h^k\)。

由于下降速度场在边界上设置为零，边界入流通量项可以忽略。

### 求解器设置

该线性系统规模预计接近

\[ 1.9\times 10^5. \]

它是标量 DG 系统，默认采用 PETSc GMRES：

- `-ksp_type gmres`
- `-pc_type bjacobi`
- `-sub_pc_type ilu`

如果单步仍然慢，可以优先降低精度，例如把 `-ksp_rtol 1e-8` 改成 `1e-7` 或 `1e-6`。

~~~freefem
string spPhi =
      "-ksp_type preonly "
    /*+ "-ksp_gmres_restart 50 "
    + "-ksp_rtol 1e-8 "
    + "-ksp_atol 1e-10 "
    + "-ksp_max_it 120 "*/
    + "-pc_type none ";
    /*+ "-sub_pc_type ilu ";*/

if (useMonitor == 1)
{
    spPhi = spPhi + "-ksp_monitor_true_residual -ksp_converged_reason ";
}

macro LFDG(phiOutput, phiInput, thetax, thetay, thetaz, T)
{
    Phih u, v;
    real dtPhi = T / levelSubSteps;
    cout << "dt = " << dtPhi << endl;
    object << "dt = " << dtPhi << endl;
    object.flush; 

    varf APhi(u, v) =
    int3d(Th)(
        (u / dtPhi + ([thetax, thetay, thetaz]' * grad(u))) * v
    )
    - intallfaces(Th)(
        (nElementonB - 1)
        * v
        * (
              0.5 * abs(N.x * thetax + N.y * thetay + N.z * thetaz)
            - 0.5 * (N.x * thetax + N.y * thetay + N.z * thetaz)
          )
        * jump(u)
    );
    varf MassPhi(u, v) = int3d(Th)( (u / dtPhi) * v );

    matrix MPhi = MassPhi(Phih, Phih);
    matrix APhim = APhi(Phih, Phih);
    Mat APhiPetsc(APhim), MPhiPetsc(MPhi);
    set(APhiPetsc, sparams = "-pc_type none");

    for (int sub = 0; sub < levelSubSteps; sub++)
    {
        v[] = MPhiPetsc * phiInput[];
        phiOutput[] = APhiPetsc^-1 * v[];
        phiInput = phiOutput;
    }
}//
~~~

---

## 9. 重初始化与迭代信息输出

### 9.1 重初始化 `RedistanceConvert`

水平集函数经过多次对流后，一般不再保持距离函数性质。这里周期性调用 FreeFEM 的 `distance` 命令，将 \(\phi\) 重新初始化为近似有向距离函数。

代码流程是：

\[ \phi_{\mathrm{DG}}\in P_1^{\mathrm{dc}} \quad\longrightarrow\quad \phi_{P1}\in P_1 \quad\longrightarrow\quad \phi_{\mathrm{Re}}\in P_1 \quad\longrightarrow\quad \phi\in P_1^{\mathrm{dc}}. \]

这里 `distance` 的作用是保持界面位置基本不变，同时恢复距离函数性质。

### 9.2 输出宏 `OutPutData`

每次迭代输出：

- 当前迭代编号；
- 固定的 `Titer`；
- 当前目标函数 `Drag`；
- 与上一步相比的目标函数变化量；
- 当前体积 `Vol`；
- 相对初始体积的变化量。

~~~freefem
macro RedistanceConvert(phiInput)
{
    P1h phiconvert, phiRe;
    phiconvert = phiInput;
    distance(Th, phiconvert, phiRe[]);
    phiInput = phiRe;
}//

macro OutPutData(i, drag, vol, obj)
{
    if (i == 0)
    {
        obj << "Initial condition" << endl;
        obj << "Drag = " << drag[i] << endl;
        obj << "Vol  = " << vol[i] << endl;
    }
    else
    {
        obj << i << "-th iteration" << endl;
        obj << "Titer = " << Titer << endl;
        obj << "Drag = " << drag[i]
            << "; diff = " << drag[i] - drag[i - 1] << endl;
        obj << "Vol  = " << vol[i]
            << "; diff from initial = " << vol[i] - vol[0] << endl;
    }

    obj.flush;
}//
~~~

---

## 10. 初始状态求解

在主迭代前，需要先完成第 0 步计算：

1. 根据初始密度 \(\rho^0=H_\gamma(\phi^0)\) 求解 Brinkman-Stokes 方程；
2. 根据 \((u^0,\rho^0,\phi^0)\) 求解下降速度场 \(\theta^0\)；
3. 计算初始目标函数

\[ J^0=J(\rho^0,u^0), \]

以及初始体积

\[ \operatorname{Vol}^0= \int_D H_\gamma(\phi^0)\,dx. \]

同时输出

\[ \operatorname{GetdJ}(\theta^0), \qquad \operatorname{GetdV}(\theta^0). \]

其中 `GetdV` 应该接近 0，用来检查体积约束是否正常。

~~~freefem
real[int] drag(iterMax + 1);
real[int] Vol(iterMax + 1);

int iter = 0;

cout << "Initial Stokes solve..." << endl;
SolveStokes(ux, uy, uz, Rho);

cout << "Initial descent velocity solve..." << endl;
SolveVelocity(thetax, thetay, thetaz, ux, uy, uz, Rho, phi);

drag[0] = GetDrag(ux, uy, uz, Rho);
Vol[0]  = GetVol(Rho);

OutPutData(0, drag, Vol, object);

real dJtheta0 = GetdJ(ux, uy, uz, thetax, thetay, thetaz, Rho, phi);
real dVtheta0 = GetdV(thetax, thetay, thetaz, phi);

object << "Initial dJ(theta) = " << dJtheta0 << endl;
object << "Initial dV(theta) = " << dVtheta0 << endl;
object.flush;

cout << "Initial drag = " << drag[0] << endl;
cout << "Initial Vol  = " << Vol[0] << endl;
~~~

---

## 11. 主迭代

主迭代从 \(n=0\) 开始，重复以下步骤：

### 第一步：水平集演化

由当前 \(\phi^n\) 和 \(\theta^n\) 解 DG 对流方程，得到

\[ \phi^{n+1}. \]

这里使用固定时间长度

\[ T_{\mathrm{iter}}, \]

不判断目标函数是否下降。

### 第二步：周期性重初始化

当迭代次数满足

\[ (n+1)\bmod \texttt{redistEvery}=0 \]

时，执行重初始化。

### 第三步：更新密度

\[ \rho^{n+1}=H_\gamma(\phi^{n+1}). \]

### 第四步：重新求解状态方程

根据新密度求

\[ (u^{n+1},p^{n+1}). \]

### 第五步：计算目标函数和体积

\[ J^{n+1}=J(\rho^{n+1},u^{n+1}), \]

\[ \operatorname{Vol}^{n+1}=\int_D \rho^{n+1}\,dx. \]

### 第六步：求下一步下降速度场

根据新的 \((u^{n+1},\rho^{n+1},\phi^{n+1})\)，求

\[ \theta^{n+1}. \]

### 第七步：按频率保存结果

每隔 `saveEvery` 步保存一次 `.vtu` 文件，用于 ParaView 后处理。

这一版没有下降性检测，所以即使某一步目标函数增大，也仍然接受更新。

~~~freefem
Phih phitest;
real timeIter = clock();

while (iter < iterMax)
{
    if (iter < 50)
        Titer = 0.0001;
    else if (50 <= iter && iter < 100)
        Titer = 0.00004;
    else if (100 <= iter && iter < 200)
        Titer = 0.000005;
    else if (200 <= iter && iter < 300)
        Titer = 0.000001;    
    cout << endl;
    cout << "iter = " << iter + 1 << endl;

    object << endl;
    object << "iter = " << iter + 1 << endl;

    phitest = phi;
    LFDG(phi, phitest, thetax, thetay, thetaz, Titer);

    if ((iter + 1) % redistEvery == 0)
    {
        RedistanceConvert(phi);
        cout << "Redistance done." << endl;
        object << "Redistance done." << endl;
    }

    Rho = HeaviSide(phi);

    SolveStokes(ux, uy, uz, Rho);

    iter++;

    drag[iter] = GetDrag(ux, uy, uz, Rho);
    Vol[iter]  = GetVol(Rho);

    OutPutData(iter, drag, Vol, object);

    SolveVelocity(thetax, thetay, thetaz, ux, uy, uz, Rho, phi);

    real dJtheta = GetdJ(ux, uy, uz, thetax, thetay, thetaz, Rho, phi);
    real dVtheta = GetdV(thetax, thetay, thetaz, phi);

    object << "dJ(theta) = " << dJtheta << endl;
    object << "dV(theta) = " << dVtheta << endl;

    cout << "Drag = " << drag[iter] << endl;
    cout << "Vol  = " << Vol[iter] << endl;
    cout << "dJ(theta) = " << dJtheta << endl;
    cout << "dV(theta) = " << dVtheta << endl;

    plot(Rho, fill = 1, value = 1, cmm = "Rho_" + iter, wait = false);

    if (iter % saveEvery == 0)
    {
        string filename = "./figure/phi_" + iter + ".vtu";
        savevtk(filename, Th, phi, Rho, [ux, uy, uz]);
    }

    object << "iteration time = " << clock() - timeIter << endl;
    cout << "iteration time = " << clock() - timeIter << endl;

    timeIter = clock();
    object.flush;
}
~~~

---

## 12. 最终输出

主循环结束后，程序保存：

1. 最终形状和流场：`FinalShape.vtu`；
2. 目标函数和体积变化数据：`Model3_Shape1.m`。

Matlab 文件中包含：

\[ \texttt{Volume(i)}, \qquad \texttt{drag(i)}. \]

可以直接在 Matlab 中绘制体积变化和目标函数变化。

~~~freefem
savevtk(ModelResultName, Th, phi, Rho, [ux, uy, uz]);

ofstream OutMatlab(ModelMatLabName);
OutMatlab << "Num = " << iter + 1 << endl;

for (int i = 0; i <= iter; i++)
{
    OutMatlab << "Volume(" << i + 1 << ") = " << Vol[i] << "; ";
    OutMatlab << "drag(" << i + 1 << ") = " << drag[i] << ";" << endl;
}

object << endl;
object << "Finished." << endl;
object << "Final drag = " << drag[iter] << endl;
object << "Final Vol  = " << Vol[iter] << endl;
object.flush;

cout << "Finished." << endl;
cout << "Final drag = " << drag[iter] << endl;
cout << "Final Vol  = " << Vol[iter] << endl;
~~~

---

## 13. 当前版本的计算量与求解器建议

你给出的规模为：

\[ \dim(S_h)=209810, \qquad \dim(\Theta_h)=27150. \]

预计 level-set DG 方程规模为

\[ \dim(\Phi_h)\approx 190000\sim 200000. \]

每次优化迭代主要耗时来自三部分：

1. Brinkman-Stokes 方程：规模约 `209810`，混合鞍点系统，最贵；
2. 下降速度场：规模约 `27150`，但需要解两次 SPD 系统；
3. Level-set DG 方程：规模约 `190000--200000`，标量 DG 对流系统。

为了尽量让每次迭代控制在 2 分钟内，本代码默认采取：

- 固定 `Titer`，不做回退试算；
- `levelSubSteps = 1`，每步只解一次 DG 方程；
- 不在每次迭代中调用 `plot`；
- `saveEvery = 10`，每 10 步保存一次；
- Stokes 用 PETSc fieldsplit；
- 下降速度场避免鞍点系统，用两个 SPD 解替代；
- level-set 用 GMRES + block Jacobi + ILU。

若实际运行仍然超过 2 分钟，优先调整：

- 将 Stokes 的 `-ksp_rtol 1e-6` 改为 `1e-5`；
- 将 level-set 的 `-ksp_rtol 1e-8` 改为 `1e-7` 或 `1e-6`；
- 将 `saveEvery = 20`；
- 确认关闭 `useMonitor = 0`；
- 若 HYPRE 不可用，将 Stokes 速度块预处理改成 GAMG。
