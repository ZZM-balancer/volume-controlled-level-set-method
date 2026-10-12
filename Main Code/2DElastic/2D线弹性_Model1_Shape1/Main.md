## 1. 必要的 load
~~~freefem
load "TetGen"
load "gmsh"
load "distance"
load "iovtk"
load "PETSc"
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
real Titer = 0.075;

int iterMax = 300; 

real mu = 8.;	    
real lambda = 1.;
real weak = 1e-4;
real gamma = 0.005;

ofstream object("./Output/Main.obj");

string ModelMatLabName = "matlab/Model1_Shape1.m";
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
|变量名|作用|
|:---:|:---:|
|Uh|位移场|
|Thetah|下降方向场|
|P1h|Lagrange线性元|
|Phih|$\phi$或$\rho$所在的有限元空间|
~~~freefem
// first label meanings: 1 = Neumann, 2 = Dirichlet, 3 = Force(Non Neumann)
// second label meanings: 1 = (V_x = 0), 2 = (V_y = 0)
border bdrightN1(t = -0.5, -0.05)  { x=1;  y=t; label=11; }; // right
border bdrightNN(t = -0.05, 0.05)  { x=1;  y=t; label=3; }; // right
border bdrightN2(t = 0.05, 0.5)    { x=1;  y=t; label=11; }; // right
border bdleftD(t = 0.5, -0.5)	   { x=-1; y=t; label=2; }; // left
border bdupN(t = 1, -1)		       { x=t;  y=0.5; label=12; };   // upper
border bdbottomN(t = -1, 1)		   { x=t;  y=-0.5; label=12; };  // lower

int n = 4;
mesh Th= buildmesh ( bdrightN1(9*n) + bdrightNN(3*n) + bdrightN2(9*n)
                     + bdleftD(20*n)
                     + bdupN(40*n)
                     + bdbottomN(40*n) );

savevtk("Th1_1.vtu", Th);

fespace Uh(Th, [P2, P2]);      
fespace Thetah(Th, [P2, P2]);   
fespace P1h(Th, P1); 
fespace Phih(Th, P1dc);

P1h h1 = hTriangle;
real hmax = h1[].max; 
real hmin = h1[].min;
object << "maximal size of a mesh triangle h = " << hmax << endl;
object << "minimal size of a mesh triangle h = " << hmin << endl;
~~~
| 网格 |
|:-----:|
|![](Th1_1.png)|

## 5. 初始形状和 PDE 求解
### 5.1 初始形状设置
先将区域分成正负值，构造func变量phi0  
再通过内置distance函数，将phi0表示的区域转换为：线性元变量phiInit.  
再将线性元变量phiInit转换为Phih变量phi  
输出密度函数$\rho = (\phi < 0) + \texttt{weak} \cdot (\phi > 0)$对应的图像
~~~freefem
macro R1()
(
    0.2 - sqrt((abs(x))^2 + (abs(y) - 0.5)^2)
) //
macro R2()
(
    0.2 - sqrt((abs(x) - 0.6)^2 + (abs(y) - 0.5)^2)
) //
macro R3()
(
    0.1 - sqrt((abs(x) - 0.3)^2 + (abs(y) - 0.25)^2)
) //
macro R4()
(
    0.1 - sqrt((abs(x) - 1)^2 + (abs(y) - 0.25)^2)
) //
macro R5()
(
    0.2 - sqrt((abs(x))^2 + (abs(y))^2)
) //
macro R6()
(
    0.2 - sqrt((abs(x) - 0.6)^2 + (abs(y))^2)
) //
func phi0 = (R1 >= 0) * R1 
            + (R2 >= 0) * R2
            + (R3 >= 0) * R3
            + (R4 >= 0) * R4
            + (R5 >= 0) * R5
            + (R6 >= 0) * R6
            + (R1 < 0 && R2 < 0 && R3 < 0 && R4 < 0 && R5 < 0 && R6 < 0) * (max(R1,max(R2,max(R3,max(R4,max(R5,R6))))));

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho;
phi = phi0;
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
|[ux, uy]|位移|
|[vx, vy]|测试函数|
|Elasticity|线弹性刚度矩阵|
|ElasticityRHS|线弹性右端项|
|ElasticityPDE|求解线弹性方程|
|Elasticityb|右端项向量|
|A|线弹性系统矩阵|


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
Uh [ux, uy], [vx, vy];

func gx = 0.;
func gy = -5.;

varf Elasticity([ux, uy], [vx, vy]) =
int2d(Th)(
        2 * mu * Rho * ( epsilon(ux, uy)' * epsilon(vx, vy) )
        + lambda * Rho * div(ux, uy) * div(vx, vy)
        )
+ on(2, ux = 0, uy = 0);

varf ElasticityRHS([ux, uy], [vx, vy]) =
int1d(Th, 3)( gx * vx + gy * vy )
+ on(2, ux = 0, uy = 0);

real[int] Elasticityb = ElasticityRHS(0, Uh);

macro ElasticityPDE(ux, uy, Rho)
{
    matrix Am = Elasticity(Uh, Uh);

    Mat A(Am);
    set(A, sparams = "-pc_type lu");

    real[int] sol = A^-1 * Elasticityb;
    ux[] = sol;
}//
~~~

### 5.3 下降方向计算
考虑 compliance 目标泛函，$\rho=H_{\gamma,s}(\phi)$：
$$
J(\Omega):= \int_D \rho A \epsilon(u):\epsilon(u) dx
$$
其中形状$\Omega = \lbrace \phi < 0 \rbrace$
对应的形状导数为：
$$
dJ(\Omega;\theta) = \int_D A\epsilon(u):\epsilon(u) \theta \cdot \nabla \rho dx
$$
在设计区域 $D$ 中，下降方向 $\theta$ 通过带体积约束的椭圆方程计算：
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
Thetah [thetax, thetay], [Wx, Wy];

varf VelocityA([thetax, thetay], [Wx, Wy]) =
int2d(Th)( 
            Grad(thetax, thetay)' * Grad(Wx, Wy)
            + ([thetax, thetay]' * [Wx, Wy])
        )
+ on(3, thetax = 0, thetay = 0)
+ on(11, 2, thetax = 0)
+ on(12, thetay = 0);
matrix A2m = VelocityA(Thetah, Thetah);

macro VelocityPDE(thetax, thetay, ux, uy, Rho)
{
    varf VelocitycRho([thetax, thetay], [Wx, Wy]) = 
    int2d(Th) ([Wx, Wy]' * grad(Rho))
    + on(3, thetax = 0, thetay = 0)
    + on(11, 2, thetax = 0)
    + on(12, thetay = 0);

    varf VelocitybRho([thetax, thetay], [Wx, Wy]) =
    - int2d(Th)( 
                 (2 * mu *  epsilon(ux, uy)' * epsilon(ux, uy) 
                   + lambda * div(ux, uy) * div(ux, uy) ) 
                 *([Wx, Wy]' * grad(Rho))
               )
    + on(3, thetax = 0, thetay = 0)
    + on(11, 2, thetax = 0)
    + on(12, thetay = 0);

    real[int] cRho = VelocitycRho(0, Thetah);
    real[int] bRho = VelocitybRho(0, Thetah);

    matrix LHSm = [[A2m, cRho], [cRho', 0]];
    real[int] RHS = [bRho, 0];
    set(LHSm, solver = CG, eps = 1e-8);
    Mat LHS(LHSm);
    set(LHS, sparams = "-pc_type lu");

    real[int] sol(RHS.n);
    sol = LHS^-1 * RHS;
    thetax[] = sol(0:sol.n-1);
}//
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
代码中下降方向在相应边界满足零位移、零分量或零法向约束，因此外边界入流项不再单独装配。
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
    Mat B(Bm), A(Am);

    for (t = 0; t < Tsum; t += dt)
    {
        v[] = B * phiInput[];
        phiOutput[] = A^-1 * v[];
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

	plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter);
    
    string filename = "./figure/Rho_"+iter+".vtu";
    string filename2 = "./figure/phi_"+iter+".vtu";
    savevtk(filename, Th, Rho);
    //savevtk(filename2, Th, phi);
    
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

