## 1. 必要的 load
本程序只使用 FreeFEM 的顺序 MUMPS 直接解法。线弹性方程、速度场方程以及 DG 更新中的线性系统均通过 `load "MUMPS_seq"` 后的 `sparsesolver` 求解，不再使用 PETSc 的 `Mat` 与 `sparams` 设置。

~~~freefem
load "TetGen"
load "gmsh"
load "distance"
load "iovtk"
load "MUMPS_seq"
verbosity = 0;
~~~

## 2. 参数设置
|变量名|作用|
|:---:|:---:|
|Titer|迭代时间步长|
|iterMax|最大迭代步数|
|mu|拉梅常数$\mu$|
|lambda|拉梅常数$\lambda$|
|gamma|界面宽度$\gamma$|
|weak|非材料区域密度$s$|
|object|输出信息保存文件(.obj)|
|ModelMatLabName|需要使用Matlab处理的数据(.m)|

~~~freefem
real Titer = 0.008;

int iterMax = 1400; 

real mu = 8.;	    
real lambda = 1.;
real weak = 1e-4;
real gamma = 0.005;

ofstream object("./Output/Main.obj");

string ModelMatLabName = "matlab/Model3_Shape1.m";
~~~

## 3. 函数和微分算子
$$
H_{\gamma,s} = (1-s)H_{\gamma}+s
$$
$$
H_{\gamma}(x) = 
\begin{cases}
	0, & \phi > \gamma, \\
    \frac{1}{2} - \frac{\phi}{\gamma} + \frac{\phi^2}{2\gamma^2}, & \phi \in (0, \gamma], \\
	\frac{1}{2} - \frac{\phi}{\gamma} - \frac{\phi^2}{2\gamma^2}, & \phi \in [-\gamma, 0], \\  
	1, & \phi < -\gamma.
\end{cases}
$$
~~~freefem
macro HeaviSide(phi)
(
    (phi > gamma) * weak + 
    (0 < phi && phi <= gamma) *
    ( (1-weak) * (0.5-phi/gamma+0.5*phi^2/gamma^2) + weak) +
    (-gamma <= phi && phi <= 0) * 
    ( (1-weak) * (0.5-phi/gamma-0.5*phi^2/gamma^2) + weak) + 
    (phi < -gamma) * 1.0
)
~~~

指标函数$sign(x)$
~~~freefem
macro sign(x)
(
    (x >= 0) * 1. + (x < 0) * (-1)
)//
~~~

$$
\begin{align*}
\epsilon(u):\epsilon(v) &= \left[\partial_x u_x, \partial_y u_y, \frac{\partial_y u_x + \partial_x u_y}{\sqrt{2}}\right] \cdot \left[\partial_x v_x, \partial_y v_y, \frac{\partial_y v_x + \partial_x v_y}{\sqrt{2}}\right] \\
\nabla u &= [\partial_x u, \partial_y u] \\
Du &= [\nabla u_x, \nabla u_y] \\
\texttt{div} \ u &= \partial_x u_x + \partial_y u_y
\end{align*}
$$

~~~freefem
macro epsilon(ux, uy) [dx(ux), dy(uy), (dy(ux)+dx(uy))/sqrt(2)] //
macro grad(u) [dx(u), dy(u)] //
macro Grad(ux, uy) [dx(ux), dy(ux), dx(uy), dy(uy)] //
macro div(ux, uy) ( dx(ux) + dy(uy) ) //
~~~

## 4. 二维网格和有限元空间
这里采用矩形设计区域 $D=[-1,1]\times[-0.5,0.5]$。为了模拟简支梁支承，底边左右两端各保留一小段支承边界：左下小线段施加 $u_x=u_y=0$，右下小线段施加 $u_y=0$。与原始文件相比，这里去掉了左侧竖边和右侧竖边上的小支承段，只保留下边界两侧的小支承段。

速度场边界仍按照几何边界分类：在受力边界 $NN$ 上令 $\theta_x=\theta_y=0$；在左右边界上令 $\theta_x=0$；在上下边界中除 $NN$ 外的部分令 $\theta_y=0$。底边两侧的小支承段也属于上下边界的一部分，因此速度场在这两段上施加 $\theta_y=0$。

|变量名|作用|
|:---:|:---:|
|Uh|位移场所在的向量 Lagrange 线性元空间 $[P_1,P_1]$|
|Thetah|速度场所在的向量 Lagrange 线性元空间 $[P_1,P_1]$|
|P1h|连续线性元空间，用于距离函数和辅助变量|
|Phih|间断线性元空间，用于 $\phi$ 和 $\rho$ 的 DG 更新|
|labNN|受力边界，同时速度场满足 $\theta_x=\theta_y=0$|
|labLR|左右边界，速度场满足 $\theta_x=0$|
|labTB|上下边界中除受力边界和支承小线段外的部分，速度场满足 $\theta_y=0$|
|labSuppL|底边左端小支承段，线弹性位移满足 $u_x=u_y=0$|
|labSuppR|底边右端小支承段，线弹性位移满足 $u_y=0$|

~~~freefem
int labNN = 3;
int labLR = 11;
int labTB = 12;
int labSuppL = 22;
int labSuppR = 222;

// labNN: force boundary NN, also theta_x = theta_y = 0
// labLR: left and right boundaries, theta_x = 0
// labTB: top and bottom boundaries except NN and support segments, theta_y = 0
// labSuppL: bottom-left support segment, ux = uy = 0 for elasticity, theta_y = 0 for velocity
// labSuppR: bottom-right support segment, uy = 0 for elasticity, theta_y = 0 for velocity
border bdbottomSuppL(t = -1, -0.95) { x = t; y = -0.5; label = labSuppL; };
border bdbottomL(t = -0.95, -0.05) { x = t; y = -0.5; label = labTB; };
border bdbottomNN(t = -0.05, 0.05) { x = t; y = -0.5; label = labNN; };
border bdbottomR(t = 0.05, 0.95)  { x = t; y = -0.5; label = labTB; };
border bdbottomSuppR(t = 0.95, 1)  { x = t; y = -0.5; label = labSuppR; };
border bdright(t = -0.5, 0.5)      { x = 1; y = t; label = labLR; };
border bdup(t = 1, -1)             { x = t; y = 0.5; label = labTB; };
border bdleft(t = 0.5, -0.5)       { x = -1; y = t; label = labLR; };

int n = 4;
mesh Th = buildmesh(
      bdbottomSuppL(n)
    + bdbottomL(18*n)
    + bdbottomNN(3*n)
    + bdbottomR(18*n)
    + bdbottomSuppR(n)
    + bdright(20*n)
    + bdup(40*n)
    + bdleft(20*n)
);

savevtk("Th2_1.vtu", Th);

fespace Uh(Th, [P1, P1]);
fespace Thetah(Th, [P1, P1]);
fespace P1h(Th, P1);
fespace Phih(Th, P1dc);

P1h h1 = hTriangle;
real hmax = h1[].max;
real hmin = h1[].min;
object << "maximal size of a mesh triangle h = " << hmax << endl;
object << "minimal size of a mesh triangle h = " << hmin << endl;

int Nu = Uh.ndof;
int Ntheta = Thetah.ndof;
~~~
| 网格 |
|:-----:|
|![](Th2_1.png)|

## 5. 初始形状和 PDE 求解
### 5.1 初始形状设置
初始形状采用两行并排圆孔，每行 $5$ 个，总计 $10$ 个孔。孔心均放在区域内部，避免与四个角点支承位置相接触。这里取圆孔半径 `holeR = 0.178`，十个孔的总体积约为 $10\pi(0.178)^2$，相对于整体区域体积 $2$，材料区域体积约为整体体积的 $50\%$。

先构造符号距离函数的初始符号函数 `phi0`，其中 $\phi>0$ 表示孔洞区域，$\phi<0$ 表示材料区域。随后通过 `distance` 重新初始化为距离函数，再转换到 DG 空间中，并输出初始密度函数 $\rho=H_{\gamma,s}(\phi)$。
~~~freefem
real holeR = 0.17;

func R11 = holeR - sqrt((x + 0.75)^2 + (y - 0.22)^2);
func R12 = holeR - sqrt((x + 0.375)^2 + (y - 0.22)^2);
func R13 = holeR - sqrt(x^2 + (y - 0.22)^2);
func R14 = holeR - sqrt((x - 0.375)^2 + (y - 0.22)^2);
func R15 = holeR - sqrt((x - 0.75)^2 + (y - 0.22)^2);

func R21 = holeR - sqrt((x + 0.75)^2 + (y + 0.22)^2);
func R22 = holeR - sqrt((x + 0.375)^2 + (y + 0.22)^2);
func R23 = holeR - sqrt(x^2 + (y + 0.22)^2);
func R24 = holeR - sqrt((x - 0.375)^2 + (y + 0.22)^2);
func R25 = holeR - sqrt((x - 0.75)^2 + (y + 0.22)^2);

func phi0 = max(
                max(max(max(max(R11, R12), R13), R14), R15),
                max(max(max(max(R21, R22), R23), R24), R25)
            );

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho;
phi = phiInit;
Rho = HeaviSide(phi);
plot(Rho, fill = 1, value = 1, cmm = "Rho_0");
savevtk("./figure/Rho_0.vtu", Th, Rho);
savevtk("./figure/phi_0.vtu", Th, phi);
~~~
| 初始形状 |
|:-----:|
|![](Rho_0.png)|

### 5.2 位移计算
|变量名|作用|
|:---:|:---:|
|[ux, uy]|位移场，使用 $[P_1,P_1]$ 元|
|[vx, vy]|测试函数，使用 $[P_1,P_1]$ 元|
|Elasticity|线弹性刚度矩阵对应的双线性型|
|ElasticityRHS|底边中部受力边界 $NN$ 上的载荷项|

在设计区域 $D$ 中，位移满足密度加权线弹性方程。受力边界为底边中部的 $NN$，底边左端小线段施加 $u_x=u_y=0$，底边右端小线段施加 $u_y=0$。这里采用小边界段近似简支点支承，但不再在左右竖边上设置额外的小支承段。

对应的变分形式为
$$
\int_D \rho A\epsilon(u):\epsilon(v)\,dx
=
\int_{NN} g\cdot v\,ds,
$$
其中 $g=(0,-20)$。线性系统使用顺序 MUMPS 直接法求解。
~~~freefem
Uh [ux, uy], [vx, vy];

func gx = 0.;
func gy = -20.;

varf Elasticity([ux, uy], [vx, vy]) =
int2d(Th)(
        2.0 * mu * Rho * (epsilon(ux, uy)' * epsilon(vx, vy))
        + lambda * Rho * div(ux, uy) * div(vx, vy)
        )
+ on(labSuppL, ux = 0, uy = 0)
+ on(labSuppR, uy = 0);

varf ElasticityRHS([ux, uy], [vx, vy]) =
int1d(Th, labNN)(gx * vx + gy * vy)
+ on(labSuppL, ux = 0, uy = 0)
+ on(labSuppR, uy = 0);

macro ElasticityPDE(ux, uy, Rho)
{
    matrix Am = Elasticity(Uh, Uh);
    real[int] Elasticityb = ElasticityRHS(0, Uh);

    set(Am, solver = sparsesolver);

    real[int] sol = Am^-1 * Elasticityb;
    ux[] = sol;
} //
~~~

### 5.3 下降方向计算
考虑目标泛函，$\rho=H_{\gamma,s}(\phi)$：
$$
J(\Omega):= \int_D \rho A\epsilon(u):\epsilon(u) dx .
$$
其中形状 $\Omega=\{\phi<0\}$。对应的形状导数写成扩展区域形式：
$$
dJ(\Omega;\theta)=\int_D A\epsilon(u):\epsilon(u)\,\theta\cdot\nabla\rho\,dx .
$$

速度场使用向量 $[P_1,P_1]$ 元离散。边界条件按照简支梁优化模型的几何边界设置：在受力边界 $NN$ 上令 $\theta_x=\theta_y=0$；在左右边界上令 $\theta_x=0$；在上下边界中除 $NN$ 之外的部分令 $\theta_y=0$。底边左右两侧的小支承段也属于下边界的一部分，因此速度场在这两段上同样取 $\theta_y=0$。这样可以使速度场在外边界上没有法向穿透，从而 DG 更新时边界通量项保持为零。

为了保持体积，速度场仍通过一个体积约束乘子修正，离散系统为
$$
\begin{bmatrix}
A_\theta & c_\rho \\
c_\rho^T & 0
\end{bmatrix}
\begin{pmatrix}
\theta \\
p
\end{pmatrix}
=
\begin{pmatrix}
b_\rho \\
0
\end{pmatrix}.
$$
该鞍点系统使用顺序 MUMPS 直接法求解。

~~~freefem
Thetah [thetax, thetay], [Wx, Wy];

varf VelocityA([thetax, thetay], [Wx, Wy]) =
int2d(Th)( 
            Grad(thetax, thetay)' * Grad(Wx, Wy)
            + ([thetax, thetay]' * [Wx, Wy])
        )
+ on(labNN, thetax = 0, thetay = 0)
+ on(labLR, thetax = 0)
+ on(labTB, labSuppL, labSuppR, thetay = 0);

matrix A2m = VelocityA(Thetah, Thetah);

macro VelocityPDE(thetax, thetay, ux, uy, Rho)
{
    varf VelocitycRho([thetax, thetay], [Wx, Wy]) = 
    int2d(Th)([Wx, Wy]' * grad(Rho))
    + on(labNN, thetax = 0, thetay = 0)
    + on(labLR, thetax = 0)
    + on(labTB, labSuppL, labSuppR, thetay = 0);

    varf VelocitybRho([thetax, thetay], [Wx, Wy]) =
    - int2d(Th)(
                 (2.0 * mu * epsilon(ux, uy)' * epsilon(ux, uy)
                   + lambda * div(ux, uy) * div(ux, uy))
                 * ([Wx, Wy]' * grad(Rho))
               )
    + on(labNN, thetax = 0, thetay = 0)
    + on(labLR, thetax = 0)
    + on(labTB, labSuppL, labSuppR, thetay = 0);

    real[int] cRho = VelocitycRho(0, Thetah);
    real[int] bRho = VelocitybRho(0, Thetah);

    matrix LHSm = [[A2m, cRho], [cRho', 0]];
    real[int] RHS = [bRho, 0];

    set(LHSm, solver = sparsesolver);

    real[int] sol = LHSm^-1 * RHS;
    thetax[] = sol(0:Ntheta-1);
} //
~~~

计算目标泛函$J(\Omega)$，形状导数$dJ(\Omega;\theta)$以及当前区域体积Vol和$d\texttt{Vol}(\Omega;\theta)$：
~~~freefem
macro GetCompliance(ux, uy, Rho)
(
    int2d(Th)(
                2 * mu * Rho * ( epsilon(ux, uy)' * epsilon(ux, uy) )
                + lambda * Rho * div(ux, uy) * div(ux, uy)
            )
) //

macro GetVol(Rho)
(
    int2d(Th) (Rho)
)//

macro GetdJ(ux, uy, thetax, thetay, Rho)
(
    int2d(Th)( 
                (2 * mu *  epsilon(ux, uy)' * epsilon(ux, uy) 
                   + lambda * div(ux, uy) * div(ux, uy) )  
                * ([thetax, thetay]' * grad(Rho))
            )
)//

macro GetdV(ux, uy, thetax, thetay, Rho)
(
    int2d(Th)( 
                [thetax, thetay]' * grad(Rho)
            )
)//
~~~

### 5.4 形状更新
采用 dual-P^1 DG方法求解：
$$
\phi_t + \theta \cdot \nabla \phi = 0
$$
在每次形状更新中取 $dt=Tsum/5$：
$$
\begin{align*}
& \int_D \frac{\phi^{n+1}-\phi^n}{dt} v dx + \int_D v \theta \cdot \nabla \phi dx - \sum_{K \in \mathcal{T}_h} \sum_{E \subset \partial K} \int_E \frac{|\theta \cdot n|-\theta \cdot n}{2}[\phi^{n+1}] v ds \\
= & \sum_{E \subset \partial D} \int_{E \bigcap \lbrace \theta \cdot n < 0 \rbrace} |\theta \cdot n | \phi^{n+1} v ds
\end{align*}
$$
由于速度场边界条件使得 $\theta\cdot n=0$ 在外边界上成立，DG 更新中的外边界通量项为零。该 DG 线性系统同样使用顺序 MUMPS 直接法求解。
~~~freefem
macro nn()
(
    N.x * thetax + N.y * thetay
)//

macro LFDG(phiOutput, phiInput, thetax, thetay, Tsum)
{
    real t = 0.;
    real dt = Tsum / 5.;

    Phih u, v, phitt, phie;
    P1h OutRho;
    
    varf Adual(u, v) =
    int2d(Th) ( (u / dt + ([thetax, thetay]' * grad(u))) * v )
    - intalledges(Th) (
                        (nTonEdge - 1) * v 
                        * (0.5 * abs(nn) - nn / 2) * jump(u) 
                    );
    
    varf DGRHS(u, v) = int2d(Th) ( (u / dt) * v );
         
    matrix Bm = DGRHS(Phih, Phih);
    matrix Am = Adual(Phih, Phih);
    set(Am, solver = sparsesolver);

    for (t = 0; t < Tsum; t += dt)
    {
        v[] = Bm * phiInput[];
        phiOutput[] = Am^-1 * v[];
        phiInput = phiOutput;
    }
}//
~~~

## 6. 程序主要部分
### 6.1 第 i 步迭代信息输出
~~~freefem
macro OutPutData(i, compliance, vol, obj)
{
    if (i == 0)
    {
        obj << "Initial Condition: " << endl;
        obj << "Initial compliance = " << compliance[i] << endl;
        obj << "Initial Vol = " << vol[i] << endl;
    }
    else
    {
        obj << i << "-th iteration: " << endl;
        obj << "Compliance = " << compliance[i] << ";   diff = " << compliance[i] - compliance[i-1] << endl;
        obj << "Vol = " << vol[i] << ";   diff = " << vol[i] - vol[0] << endl;
    }
    obj.flush;
}//
~~~
### 6.2 重新初始化：利用 distance 函数
~~~freefem
macro RedistanceConvert(phiInput)
{
    P1h phiRe;
    distance(Th, phiInput, phiRe[]);
    phiInput = phiRe;
}//
~~~
### 6.3 第 0 步：先计算位移和下降方向

~~~freefem
Uh [uxnew, uynew];
Phih phitest, phinew, Rhonew;
real[int] compliance(iterMax + 1);
real[int] Vol(iterMax + 1);
real compliancenew, Volnew;

int iter = 0;
cout << "Yes1!" << endl; 
ElasticityPDE(ux, uy, Rho);  
cout << "Yes2!" << endl; 
VelocityPDE(thetax, thetay, ux, uy, Rho);

compliance[iter] = GetCompliance(ux, uy, Rho);
Vol[iter] = GetVol(Rho);
OutPutData(iter, compliance, Vol, object);
~~~
### 6.4 迭代求解
~~~freefem
while (iter < iterMax)
{
    cout << "Iter = " << iter << endl;
    object << endl;
    object << "Iter = " << iter << endl;

    phitest = phi;
    LFDG(phi, phitest, thetax, thetay, Titer);
    if ((iter - 1) % 8 == 0)
    {
        RedistanceConvert(phi);
        cout << "Redistance End!" << endl;
    }
    Rho = HeaviSide(phi);
    ElasticityPDE(ux, uy, Rho);
    iter++;
    compliance[iter] = GetCompliance(ux, uy, Rho);
    Vol[iter] = GetVol(Rho); 
            
    object << "Titer = " << Titer << endl;

    OutPutData(iter, compliance, Vol, object);

	plot(Rho, fill = 1, value = 0, cmm = "Rho_"+iter);
    
    string filename = "./figure/Rho_"+iter+".vtu";
    string filename2 = "./figure/phi_"+iter+".vtu";
    if (iter % 10 == 0)
    {
        savevtk(filename, Th, Rho);
        savevtk(filename2, Th, phi);
    }
    
    VelocityPDE(thetax, thetay, ux, uy, Rho);
}
~~~

### 6.5 compliance 信息和体积信息输出（输出为 Matlab 文件）
~~~freefem
ofstream OutMatlab(ModelMatLabName);
OutMatlab << "Num = " << iter + 1 << endl;
for (int i = 0; i <= iter; i++)
{
    OutMatlab << "Volume(" << i+1 << ") = " << Vol[i] << "; ";
    OutMatlab << "compliance(" << i+1 << ") = " << compliance[i] << ";"<<endl;
}
~~~

## 7. 实验结果
### 7.1 中间形状展示
|| 形状 |
|:----:|:-----:|
|$\rho_0$|![](Rho_0.png)|
|$\rho_{20}$|![](Rho_20.png)|
|$\rho_{60}$|![](Rho_60.png)|
|$\rho_{120}$|![](Rho_120.png)|
|$\rho_{200}$|![](Rho_200.png)|
|$\rho_{\texttt{final}}$|![](Rho_300.png)|

### 7.2 compliance 变化和体积变化图像
||图像|
|:----:|:-----:|
|compliance图像|![](Compliance.png)|
|体积图像|![](Volume.png)|

