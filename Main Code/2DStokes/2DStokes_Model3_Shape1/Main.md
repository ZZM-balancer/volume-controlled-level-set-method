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
|alphaLow|$\underline{\alpha}$|
|alphaUp|$\overline{\alpha}$|
|qq|$q$|
|gamma|界面厚度$\gamma$|
|object|输出信息保存文件(.obj)|
|ModelMatLabName|需要使用Matlab处理的数据(.m)|

~~~freefem
real Titer = 0.005;

int iterMax = 300; 

real alphaLow = 2.5*1e-4;
real alphaUp  = 2.5*1e+4;
real qq = 0.1;
real gamma = 0.005;

ofstream object("./Output/Main.obj");

string ModelMatLabName = "matlab/Model3_Shape1.m";
~~~

## 3. 函数和微分算子
Heaviside 函数
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
    (phi > gamma) * 0.0 + 
    (0 < phi && phi <= gamma) *
    ( 0.5-phi/gamma+0.5*phi^2/gamma^2 ) +
    (-gamma <= phi && phi <= 0) * 
    ( 0.5-phi/gamma-0.5*phi^2/gamma^2 ) + 
    (phi < -gamma) * 1.0
)//
~~~

$$
\begin{align*}
\alpha(\rho) &= \overline{\alpha} + (\underline{\alpha} - \overline{\alpha}) \frac{\rho(1+q)}{\rho + q} \\
\alpha'(\rho) &=  (\underline{\alpha} - \overline{\alpha}) \frac{q(1+q)}{(\rho+q)^2}
\end{align*}
$$

~~~freefem
macro alpha(Rho)
(
    alphaUp + (alphaLow - alphaUp) * (qq + 1) * Rho / (qq + Rho)
)//

macro dalpha(Rho)
(
    (alphaLow - alphaUp) * (qq + 1) * qq / ((qq + Rho)*(qq + Rho))
)//
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
\nabla u &= [\partial_x u, \partial_y u] \\
Du &= [\nabla u_x, \nabla u_y] \\
\texttt{div} \ u &= \partial_x u_x + \partial_y u_y
\end{align*}
$$

~~~freefem
macro grad(u) [dx(u), dy(u)] //
macro Grad(ux, uy) [dx(ux), dy(ux), dx(uy), dy(uy)] //
macro div(ux, uy) ( dx(ux) + dy(uy) ) //
~~~

## 4. 二维网格和有限元空间
|变量名|作用|
|:---:|:---:|
|Uh|速度场|
|Thetah|下降方向场|
|P1h|Lagrange线性元|
|Phih|$\phi$或$\rho$所在的有限元空间|
|Sth|Uh $\times$ P1h|
~~~freefem
// first label meanings:  1 = Dirichlet, 2 = in, 3 = out
border bdrightwall1(t = 0, 1)         {x = 1; y = t; label = 1;}
border bdupwall(t = 1, 0)             {x = t; y = 1; label = 1;}
border bdleftwall1(t = 1, 0.9)        {x = 0; y = t; label = 1;}
border bdin(t = 0.9, 0.7)             {x = 0; y = t; label = 2;}
border bdleftwall2(t = 0.7, 0)        {x = 0; y = t; label = 1;}
border bdbottomwall2(t = 0, 0.7)      {x = t; y = 0; label = 1;}
border bdout(t = 0.7, 0.9)            {x = t; y = 0; label = 3;}
border bdbottomwall1(t = 0.9, 1)      {x = t; y = 0; label = 1;}

int n = 14;
mesh Th = buildmesh (bdrightwall1(10 * n) 
                     + bdupwall(10 * n) 
                     + bdleftwall1(n) + bdin(2 * n) + bdleftwall2(7 * n)
                     + bdbottomwall1(n) + bdout(2 * n) + bdbottomwall2(7 * n));

savevtk("Th3_1.vtu", Th);

fespace Uh(Th, [P2, P2]);      
fespace Thetah(Th, [P2, P2]);   
fespace P1h(Th, P1);           
fespace Sth = Uh * P1h;
fespace Phih(Th, P1dc);

P1h h1 = hTriangle;
real hmax = h1[].max; 
real hmin = h1[].min;
object << "maximal size of a mesh triangle h = " << hmax << endl;
object << "minimal size of a mesh triangle h = " << hmin << endl;
~~~
| 网格 |
|:-----:|
|![](Th2_1.png)|

## 5. 初始形状和 PDE 求解
### 5.1 初始形状设置
先将区域分成正负值，构造func变量phi0  
再通过内置distance函数，将phi0表示的区域转换为：线性元变量phiInit.  
再将线性元变量phiInit转换为Phih变量phi  
输出密度函数$\rho = (\phi < 0)$对应的图像
~~~freefem
macro R()
(
   max(x, y)
) //

func phi0 = (R <= 0.7) * (0.7 - R)
            + (R > 0.7 && R < 0.9) * (-min(R - 0.7, 0.9 - R))
            + (R >= 0.9) * (R - 0.9);

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho;
phi = phi0;
Rho = HeaviSide(phi);
plot(Rho, fill = 1, value = 1, cmm = "Rho_0");
savevtk("./figure/Rho_0.vtu", Th, Rho);
~~~
| 初始形状 |
|:-----:|
|![](Rho_0.png)|

### 5.2 速度-压力计算
|变量名|作用|
|:---:|:---:|
|[ux, uy]|速度|
|[vx, vy]|测试函数|
|StokesSteady|状态方程中与 $\rho$ 无关的矩阵部分|
|StokesRho|Brinkman 惩罚项 $\alpha(\rho)u\cdot v$ 对应的矩阵|
|StokesRHS|状态方程右端项及 Dirichlet 数据|
|A1m|状态方程固定矩阵部分|
|Stokesb|右端项向量|


在设计区域 $D$ 中，速度 $u$ 与压力 $p$ 满足 Brinkman-Stokes 变分形式：
$$
\begin{align*}
\begin{cases}
\int_D \nabla u : \nabla v dx - \int_D p \ \texttt{div} \ v dx + \int_D \alpha(\rho) u \cdot v dx = \int_D f \cdot v dx & \forall v \in [H_{\partial D_D}^1(D)]^2 \\
\int_D q \ \texttt{div} \ u dx = 0 & \forall q \in L_0^2(D)
\end{cases}
\end{align*}
$$
对应线性系统：
$$
(A + A_{\rho})x = b
$$
$A$ 是 Stokes 与压力耦合的固定矩阵，$A_{\rho}$ 是 $\alpha(\rho)u\cdot v$ 对应的 Brinkman 惩罚矩阵，$b$ 是体力和边界速度数据形成的右端项。
~~~freefem
P1h p, q;
Uh [ux, uy], [vx, vy];

func gx1 = 1 - (10*y - 8)^2;
func gx2 = 0;
func gy1 = 0;
func gy2 = -1 + (10*x - 8)^2;
func fx = 0;
func fy = 0;

varf StokesSteady(<[ux, uy], [p]>, <[vx, vy], [q]>) = 
int2d(Th)(
            Grad(ux, uy)' * Grad(vx, vy)
            - p * div(vx, vy)
            - q * div(ux, uy) 
            - p * q * 1e-8
        )
+ on(2, ux = gx1, uy = gy1)
+ on(1, ux = 0, uy = 0);

varf StokesRHS(<[ux, uy], [p]>, <[vx, vy], [q]>) = 
int2d(Th) (fx * vx + fy * vy)
+ on(2, ux = gx1, uy = gy1)
+ on(1, ux = 0, uy = 0);

matrix A1m = StokesSteady(Sth, Sth);
real[int] Stokesb = StokesRHS(0, Sth);

macro StokesPDE(ux, uy, Rho)
{
    matrix Am;

    varf StokesRho(<[ux, uy], [p]>, <[vx, vy], [q]>) = 
    int2d(Th)( alpha(Rho) * ([ux, uy]' * [vx, vy]) );
    matrix ARho = StokesRho(Sth, Sth);

    Am = A1m + ARho;

    Mat A(Am);
    set(A, sparams = "-pc_type lu");

    real[int] sol = A^-1 * Stokesb;
    [ux[], p[]] = sol;
}//
~~~

### 5.3 下降方向计算
考虑阻力型目标泛函，$\rho=H_\gamma(\phi)$：
$$
J(\Omega) := \frac{1}{2} \int_D \alpha(\rho) |u|^2 dx + \frac{1}{2} \int_D |\nabla u|^2 dx - \int_D f \cdot u dx
$$
其中形状$\Omega = \lbrace \phi < 0 \rbrace$
对应的形状导数为：
$$
dJ(\Omega;\theta)= \int_D -\frac{1}{2}\alpha'(\rho)|u|^2\,\theta\cdot\nabla\rho\,dx
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
+ on(1, 2, 3, thetax = 0, thetay = 0);

matrix A2m = VelocityA(Thetah, Thetah);

macro VelocityPDE(thetax, thetay, ux, uy, Rho)
{
    varf VelocitycRho([thetax, thetay], [Wx, Wy]) = 
    int2d(Th) ([Wx, Wy]' * grad(Rho))
    + on(1, 2, 3, thetax = 0, thetay = 0);

    varf VelocitybRho([thetax, thetay], [Wx, Wy]) =
    int2d(Th)( 
                dalpha(Rho) * (0.5 * ([ux, uy]' * [ux, uy]))  
                * ([Wx, Wy]' * grad(Rho))
            )
    + on(1, 2, 3, thetax = 0, thetay = 0);

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
macro GetDrag(ux, uy, Rho)
(
    int2d(Th) ( 0.5 * alpha(Rho) * ([ux, uy]' * [ux, uy]))
    + int2d(Th) ( 0.5 * (Grad(ux, uy)' * Grad(ux, uy)) )
    - int2d(Th) (fx * ux + fy * uy)
) //

macro GetVol(Rho)
(
    int2d(Th) (Rho)
)//

macro GetdJ(ux, uy, thetax, thetay, Rho)
(
    int2d(Th)( 
                ( -dalpha(Rho) * (0.5 * ([ux, uy]' * [ux, uy])) ) 
                * [thetax, thetay]' * grad(Rho)
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
代码中下降方向在相应边界满足零法向或零分量约束，因此外边界入流项不再单独装配。
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
macro OutPutData(i, drag, vol, obj)
{
    if (i == 0)
    {
        obj << "Initial Condition: " << endl;
        obj << "Initial drag = " << drag[i] << endl;
        obj << "Initial Vol = " << vol[i] << endl;
    }
    else
    {
        obj << i << "-th iteration: " << endl;
        obj << "Drag = " << drag[i] << ";   diff = " << drag[i] - drag[i-1] << endl;
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
### 6.3 第 0 步：先计算速度-压力和下降方向

~~~freefem
Uh [uxnew, uynew];
Phih phitest, phinew, Rhonew;
real[int] drag(iterMax + 1);
real[int] Vol(iterMax + 1);
real dragnew, Volnew;

int iter = 0;

StokesPDE(ux, uy, Rho);   
VelocityPDE(thetax, thetay, ux, uy, Rho);

drag[iter] = GetDrag(ux, uy, Rho);
Vol[iter] = GetVol(Rho);
OutPutData(iter, drag, Vol, object);
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
    if ((iter - 1) % 20 == 0)
    {
        RedistanceConvert(phi);
        cout << "Redistance End!" << endl;
    }
    Rho = HeaviSide(phi);
    StokesPDE(ux, uy, Rho);
    iter++;
    drag[iter] = GetDrag(ux, uy, Rho);
    Vol[iter] = GetVol(Rho); 
            
    object << "Titer = " << Titer << endl;

    OutPutData(iter, drag, Vol, object);

	plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter);
    
    string filename = "./figure/Rho_"+iter+".vtu";
    if (iter % 10 == 0)
    {
        savevtk(filename, Th, Rho, [ux, uy]);
    }
    
    VelocityPDE(thetax, thetay, ux, uy, Rho);
}
~~~

### 6.5 drag 信息和体积信息输出（输出为 Matlab 文件）
~~~freefem
ofstream OutMatlab(ModelMatLabName);
OutMatlab << "Num = " << iter + 1 << endl;
for (int i = 0; i <= iter; i++)
{
    OutMatlab << "Volume(" << i+1 << ") = " << Vol[i] << "; ";
    OutMatlab << "drag(" << i+1 << ") = " << drag[i] << ";"<<endl;
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

### 7.2 drag 变化和体积变化图像
||图像|
|:----:|:-----:|
|drag 图像|![](drag.png)|
|体积图像|![](Volume.png)|

