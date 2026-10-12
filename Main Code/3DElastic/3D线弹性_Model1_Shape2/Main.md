## 1. 必要的 load
~~~freefem
// NBPROC 4
// PRARM -ns -wg

load "TetGen"
load "gmsh"
load "distance"
load "iovtk"
load "PETSc"
/*include "macro_ddm.idp"*/
verbosity = 0;
~~~

## 2. 参数设置
|变量名|作用|
|:---:|:---:|
|Titer|单次迭代时间步长|
|Tmax|时间步长上限|
|iterMax|最大迭代步数|
|mu, lambda|拉梅常数$\mu,\lambda$|
|weak|非材料区域密度$s$|
|gamma|界面厚度$\gamma$|
|object|输出信息保存文件,.obj|
|ModelShapeName|初始形状,.vtu|
|ModelResultName|优化结果形状,.vtu|
|ModelMatLabName|需要使用Matlab处理的数据,.m|

~~~freefem
real Tmax = 0.0025;
real Titer; 
real iterMax = 300;

real mu = 10.0/26.0;	    
real lambda = 15.0/26.0;
real weak = 5e-3;
real gamma = 0.03;

ofstream object("./Output/Main.obj");

string ModelMatLabName = "matlab/Model1_Shape2.m";
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
\nabla u &= [\partial_x u, \partial_y u, \partial_z u] \\
Du &= [\nabla u_x, \nabla u_y, \nabla u_z] \\
\texttt{div} \ u &= \partial_x u_x + \partial_y u_y + \partial_z u_z
\end{align*}
$$

~~~freefem
macro epsilon(ux, uy, uz) [dx(ux), dy(uy), dz(uz), (dy(ux)+dx(uy))/sqrt(2.), (dz(ux)+dx(uz))/sqrt(2.), (dz(uy)+dy(uz))/sqrt(2.)] //
macro grad(u) [dx(u), dy(u), dz(u)] //
macro Grad(ux, uy, uz) [dx(ux), dy(ux), dz(ux), dx(uy), dy(uy), dz(uy), dx(uz), dy(uz), dz(uz)] //
macro div(ux, uy, uz) ( dx(ux) + dy(uy) + dz(uz) ) //
~~~

计算目标泛函$J(\Omega)$，(Test是右端项计算结果)

形状导数$dJ(\Omega;\theta)$

以及当前区域体积Vol

和$d\texttt{Vol}(\Omega;\theta)$：
~~~freefem
macro GetTestvalue(phi)
(
    sqrt(
    int3d(Th)(
        (sqrt( (dx(phi))^2 + (dy(phi))^2 + (dz(phi))^2 ) - 1.0)^2
    )
    )
)//

macro GetCompliance(Aeueu, Rho)
(
    int3d(Th)(
                Aeueu * Rho
            )
) //

macro GetVol(Rho)
(
    int3d(Th) (Rho)
)//

macro GetdJ(Aeueu, thetax, thetay, thetaz, Rho)
(
    int3d(Th)( 
                Aeueu * ([thetax, thetay, thetaz]' * grad(Rho))
            )
)//

macro GetdV(thetax, thetay, thetaz, Rho)
(
    int3d(Th)( 
                [thetax, thetay, thetaz]' * grad(Rho)
            )
)//
~~~

## 4. 三维网格和有限元空间
|变量名|作用|
|:---:|:---:|
|Uh|位移场|
|Thetah|下降方向场|
|P1h|Lagrange线性元|
|Phih|$\phi$或$\rho$所在的有限元空间|
|Sth|Uh $\times$ P1h|
|Veh|Thetah $\times$ P1h|

~~~freefem
mesh3 Th = gmshload3("Th006_010.msh");
/*mesh3 Th = gmshload3("Th.msh");*/
real r = 0.1;
macro ForceFace()
(
    (sqrt((x - 0.5)^2 + (z - 0.5)^2) <= r) && 
    (abs(y-2) <= 1e-6) 
)//

macro DirichletFace()
(
    (abs(y) <= 1e-6)
)//

macro NeumannFacex()
(
    (abs(x) <= 1e-6) || (abs(x-1) <= 1e-6)
)//

macro NeumannFacey()
(
    (sqrt((x - 0.5)^2 + (z - 0.5)^2) > r) && 
    (abs(y-2) <= 1e-6)
)//

macro NeumannFacez()
(
    (abs(z) <= 1e-6) || (abs(z-1) <= 1e-6)
)//

/* First label: 1--Neumann, 2--Dirichlet, 3--Neumann with force
   Second label: 1--(V_x = 0), 2--(V_y = 0), 3--(V_z = 0) */
macro Labelfunc()
(
    3 * ForceFace +
    2 * DirichletFace +
    11 * NeumannFacex +
    12 * NeumannFacey +
    13 * NeumannFacez +
    0
)// 

Th = change(Th, flabel = Labelfunc);

savevtk("ThModel1_2.vtu", Th);

fespace P0h(Th, P0);
fespace Uh(Th, [P1, P1, P1]);       
fespace Thetah(Th, [P1, P1, P1]);   
fespace P1h(Th, P1);
fespace PTh(Th, P1);           
fespace Phih(Th, P1dc3d);       
fespace Sth = Uh * P1h;

P1h h1 = hTriangle;
real hmax = h1[].max; 
real hmin = h1[].min;
object << "Number of the vertices = " << Th.nv << endl;
object << "Number of the Elements = " << Th.nt << endl;
object << "maximal size of a mesh triangle h = " << hmax << endl;
object << "minimal size of a mesh triangle h = " << hmin << endl;
object << "Number of Matrix size1 = " << Uh.ndof << endl;
object << "Number of Matrix size2 = " << Thetah.ndof << endl;
object << "Number od Matrix size3 = " << Phih.ndof << endl;
object.flush;
~~~
| 网格 |
|:-----:|
|![](Th1_1.png)|

## 5. 初始形状和 PDE 求解
### 5.1 初始形状设置
先将区域分成正负值，构造func变量phi0  
再通过内置distance函数，将phi0表示的区域转换为：线性元变量phiInit.  
再将线性元变量phiInit转换为Phih变量phi  
输出密度函数$\rho = (\phi < 0)$对应的图像
~~~freefem
real rball = 0.22;
macro Sp() (sqrt((x - 0.5)^2 + (z - 0.5)^2) - rball) //

func phi0 = 
(Sp >= 0) * Sp +
(Sp < 0) * 
(
max(Sp, -min(y, 2-y))
);

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho;
/* P1h Rho; */
phi = phiInit;
Rho = HeaviSide(phi);
plot(Rho, fill = 1, value = 1, cmm = "Rho_0");
savevtk("./figure/phi_0.vtu", Th, phi);
~~~
| 初始形状 |
|:-----:|
|![](Rho_0.png)|

### 5.2 位移计算
|变量名|作用|
|:---:|:---:|
|[ux, uy, uz]|位移${\bold{u}}=(u_x,u_y,u_z)$|
|[vx, vy, vz]|测试函数${\bold{v}}=(v_x,v_y,v_z)$|
|(gx, gy, gz)|作用力$g = (g_x,g_y,g_z)$|
|Elasticity|计算线弹性方程刚度矩阵|
|ElasticityRHS|计算线弹性方程右端项|
|A|$A$|
|Elasticityb|$b$|
|ElasticityPDE|线弹性方程求解|


在设计区域 $D$ 中，位移 $u$ 满足线弹性变分形式：
$$
\begin{align*}
\int_D \rho_{\texttt{weak}} A \varepsilon(u):\varepsilon(v) dx = \int_D f \cdot u dx + \int_{\partial D_N} g \cdot u ds, \quad \forall v \in [H_{\partial D_D}^1(D)]^2
\end{align*}
$$
对应线性系统：
$$
Ax = b
$$
$A$ 是由 $2\mu\rho\,\varepsilon(u):\varepsilon(v)+\lambda\rho\,\operatorname{div}u\operatorname{div}v$ 形成的线弹性刚度矩阵，$b$ 是边界力右端项向量。
~~~freefem
P1h p, q;
Uh [ux, uy, uz], [vx, vy, vz];

func gx = 0;
func gy = 0;
func gz = -1.0/(pi*r*r);

macro GetComplianceRHS(ux, uy, uz)
(
    int2d(Th, 3) (gx * ux + gy * uy + gz * uz)
)//

varf Elasticity([ux, uy, uz], [vx, vy, vz]) = 
int3d(Th)(
            2 * mu * Rho * ( epsilon(ux, uy, uz)' * epsilon(vx, vy, vz) )
            + lambda * Rho * div(ux, uy, uz) * div(vx, vy, vz)
        )
+ on(2, ux = 0, uy = 0, uz = 0);

varf ElasticityRHS([ux, uy, uz], [vx, vy, vz]) = 
int2d(Th, 3) (gx * vx + gy * vy + gz * vz)
+ on(2, ux = 0, uy = 0, uz = 0);

real[int] Elasticityb = ElasticityRHS(0, Uh);

macro ElasticityPDE(ux, uy, uz, Rho)
{
    real time0 = clock();
    matrix Am = Elasticity(Uh, Uh);
    Mat A(Am);
    cout << "Elasticity construct time = " << clock() - time0 << endl;

    time0 = clock();
    set(A, sparams = "-ksp_type cg -ksp_rtol 1e-10");
    real[int] sol = A^-1 * Elasticityb;
    cout << "Elasticity solve time by KSP = " << clock() - time0 << endl;

    ux[] = sol;
    
}//
~~~

### 5.3 下降方向计算
|变量名|作用|
|:---:|:---:|
|[thetax, thetay, thetaz]|速度场${\bf{\theta}}=(\theta_x,\theta_y,\theta_z)$|
|[vx, vy, vz]|测试函数|
|VelocityA|计算速度场方程刚度矩阵的$A$|
|A|$A$|
|VelocitycRho|计算速度场方程刚度矩阵的$c_{\rho}$|
|VelocitybRho|计算速度场方程右端项的$b_{\rho}$|
|ElasticityPDE|线弹性方程求解|
|LHS|速度场方程刚度矩阵|
|RHS|速度场方程右端项|

考虑 compliance 目标泛函，$\rho=H_{\gamma,s}(\phi)$：
$$
J(\Omega):= \int_D \rho A \epsilon(u):\epsilon(u) dx
$$
其中形状$\Omega = \lbrace \phi < 0 \rbrace$
对应的形状导数为：
$$
dJ(\Omega;\theta) = \int_D A\epsilon(u):\epsilon(u) \theta \cdot \nabla \rho dx
$$
在设计区域 $D$ 中，位移 $u$ 满足线弹性变分形式：
$$
\begin{align*}
\begin{cases}
\int_D \nabla \theta : \nabla w dx + \int_D  \theta \cdot w dx + \int_D p \ w \cdot \nabla \rho dx  = -dJ(\Omega;w) & \forall w \in [H_0^1(D)]^2 \\
\int_D q \ \theta \cdot \nabla \rho dx = 0 & \forall q \in \mathbb{R}
\end{cases}
\end{align*}
$$
对应的线性系统：
$$
\begin{bmatrix}
A & c_{\rho} \\
c_{\rho}^T & 0
\end{bmatrix} 
\begin{pmatrix}
V_{\rho} \\
p
\end{pmatrix} =
\begin{pmatrix}
b_{\rho} \\
0
\end{pmatrix} 
$$

~~~freefem
Thetah [thetax, thetay, thetaz], [Wx, Wy, Wz];
PTh thetaM;

varf VelocityA([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
int3d(Th)( 
            Grad(thetax, thetay, thetaz)' * Grad(Wx, Wy, Wz)
            + ([thetax, thetay, thetaz]' * [Wx, Wy, Wz])
        )
+ on(3, thetax = 0, thetay = 0, thetaz = 0)
+ on(11, thetax = 0)
+ on(2, 12, thetay = 0)
+ on(13, thetaz = 0);

matrix A2m = VelocityA(Thetah, Thetah);

macro VelocityPDE(thetax, thetay, thetaz, Aeueu, Rho)
{
    varf VelocitycRho([thetax, thetay, thetaz], [Wx, Wy, Wz]) = 
    int3d(Th) ([Wx, Wy, Wz]' * grad(Rho))
    + on(3, thetax = 0, thetay = 0, thetaz = 0)
    + on(11, thetax = 0)
    + on(2, 12, thetay = 0)
    + on(13, thetaz = 0);

    varf VelocitybRho([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
    - int3d(Th)( 
                 Aeueu * ([Wx, Wy, Wz]' * grad(Rho))
               )
    + on(3, thetax = 0, thetay = 0, thetaz = 0)
    + on(11, thetax = 0)
    + on(2, 12, thetay = 0)
    + on(13, thetaz = 0);

    real time0 = clock();   
    real[int] bRho = VelocitybRho(0, Thetah); 
    real[int] cRho = VelocitycRho(0, Thetah);
    matrix LHSm = [[A2m, cRho], [cRho', 1e-10]];    
    Mat LHS = LHSm;
    real[int] RHS = [bRho, 0];
    cout << "Velocity construct time = " << clock() - time0 << endl; 
    
    time0 = clock();   
    real[int] sol(RHS.n);
    real p;
    set(LHS, sparams = "-ksp_type cg -ksp_rtol 1e-8");
    sol = LHS^-1 * RHS;
    [thetax[], p] = sol;
    cout << "Velocity solve time by KSP = " << clock() - time0 << endl;
    cout <<  " p = " << p << endl;
    /*
    cout << thetax[].n << " " << sol.n << endl;
    cout << "dVol = " << GetdV(thetax, thetay, thetaz, Rho) << " p = " << sol(sol.n-1) << endl;
    */
}//
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
代码中下降方向在相应边界满足零法向或零分量约束，因此外边界入流项不再单独装配。
~~~freefem
macro nn()
(
    N.x * thetax + N.y * thetay + N.z * thetaz
)//

macro LFDG(phiOutput, phiInput, thetax, thetay, thetaz, Tsum)
{
    real t = 0.;
    real dt = Tsum / 5.;

    Phih u, v, Rhotest;
    P1h OutRho;
    
    varf Adual(u, v) =
    int3d(Th) ( (u / dt + ([thetax, thetay, thetaz]' * grad(u))) * v )
    - intallfaces(Th) (
                        (nElementonB - 1) * v 
                        * (0.5 * abs(nn) - nn / 2) * jump(u) 
                    );
    
    varf DGRHS(u, v) = int3d(Th) ( (u / dt) * v );
    
    real time0 = clock();
    matrix Bm = DGRHS(Phih, Phih);
    matrix Am = Adual(Phih, Phih);
    Mat B(Bm), A(Am);
    set(A, sparams = "-ksp_type gmres -ksp_rtol 1e-10");
    cout << "DG construct time = " << clock() - time0 << endl;

    for (t = 0; t < Tsum; t += dt)
    {
        time0 = clock();
        v[] = B * phiInput[];
        phiOutput[] = A^-1 * v[];
        phiInput = phiOutput;
        cout << "DG solve time = " << clock() - time0 << endl;
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
        obj << "Initial Compliance = " << compliance[i] << endl;
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
Uh [uxnew, uynew, uznew];
Phih phitest, phinew, Rhonew;
P0h Aeueu;
real[int] compliance(iterMax + 1);
real[int] Vol(iterMax + 1);
real compliancenew, Volnew;

int iter = 0;

cout << "Time Begin "  << endl;
real time = clock();
ElasticityPDE(ux, uy, uz, Rho);
Aeueu = (2 * mu *  epsilon(ux, uy, uz)' * epsilon(ux, uy, uz) 
                   + lambda * div(ux, uy, uz) * div(ux, uy, uz) );
cout << "Elasticity time = " << clock() - time << endl;

time = clock();
VelocityPDE(thetax, thetay, thetaz, Aeueu, Rho);
cout << "Velocity time = " << clock() - time << endl;
savevtk("test_V.vtu", Th, [thetax, thetay, thetaz]);

cout << endl;
object << endl;
compliance[iter] = GetComplianceRHS(ux, uy, uz);
Vol[iter] = GetVol(Rho);
OutPutData(iter, compliance, Vol, object);
~~~
### 6.4 迭代求解
~~~freefem
time = clock();
int Redis = 0;
while ( iter < iterMax )
{
    thetaM = (thetax * thetax + thetay * thetay + thetaz * thetaz)^(0.5);
    Titer = min(Tmax, 0.4*hmin/thetaM[].max);
    /*Titer = Tmax;*/

    iter++;
    Redis++;
    cout << endl;
    cout << "iter = " << iter << endl;
    object << endl;
    object << "iter = " << iter << endl;
    object << "Titer = " << Titer << endl;

    phitest = phi;
    real time0 = clock();
    LFDG(phi, phitest, thetax, thetay, thetaz, Titer);
    
    real tv = GetTestvalue(phi);
    cout << "|| |grad(phi)| -1 || = " << tv << endl;
    if (Redis >= 12 && tv >= 0.2)
    {
        RedistanceConvert(phi);
        cout << "Redistance End!" << endl;
        cout << "Redistance || |grad(phi)| -1 || = " << GetTestvalue(phi) << endl;
        Redis = 0;
    }
    Rho = HeaviSide(phi);
    cout << "LFDG time = " << clock() - time0 << endl;

    time0 = clock();
    ElasticityPDE(ux, uy, uz, Rho);
    Aeueu = (2 * mu *  epsilon(ux, uy, uz)' * epsilon(ux, uy, uz) 
                   + lambda * div(ux, uy, uz) * div(ux, uy, uz) );
    cout << "Elasticity time = " << clock() - time0 << endl;

    time0 = clock();
    VelocityPDE(thetax, thetay, thetaz, Aeueu, Rho);
    cout << "Velocity time = " << clock() - time0 << endl;
    cout << "Test boundary value = " << int2d(Th)(thetax * N.x + thetay * N.y + thetaz * N.z) << endl;
    
    object << "Test Compliance = " << GetCompliance(Aeueu, Rho) << endl;
    compliance[iter] = GetComplianceRHS(ux, uy, uz);
    Vol[iter] = GetVol(Rho); 
   
    OutPutData(iter, compliance, Vol, object);
    
    plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter);
    /*if (iter % 10 == 0)*/
    {
        string filename = "./figure/phi_"+iter+".vtu";
        savevtk(filename, Th, phi);
        filename = "./figure/Rho_"+iter+".vtu";
        savevtk(filename, Th, Rho);
    }
    
    object << "iter Time = " << clock() - time << endl;   
    object.flush;
    cout << "iter Time = " << clock() - time << endl;
    time = clock();
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
|$\rho_5$|![](Rho_5.png)|
|$\rho_{10}$|![](Rho_10.png)|
|$\rho_{15}$|![](Rho_15.png)|
|$\rho_{20}$|![](Rho_20.png)|
|$\rho_{30}$|![](Rho_30.png)|
|$\rho_{50}$|![](Rho_50.png)|
|$\rho_{\texttt{final}}$|![](Rho_final.png)|

### 7.2 compliance 变化和体积变化图像
||图像|
|:----:|:-----:|
|compliance 图像|![](drag.png)|
|体积图像|![](Volume.png)|