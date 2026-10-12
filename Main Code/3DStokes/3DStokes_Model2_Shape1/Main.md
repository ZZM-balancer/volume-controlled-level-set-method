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
|Titer|起始迭代时间步长|
|dTmax|最大迭代时间步长|
|iterMax|最大迭代步数|
|LoopMax|最大循环次数|
|alphaLow|$\underline{\alpha}$|
|alphaUp|$\overline{\alpha}$|
|qq|$q$|
|object|输出信息保存文件,.obj|
|ModelShapeName|初始形状,.vtu|
|ModelResultName|优化结果形状,.vtu|
|ModelMatLabName|需要使用Matlab处理的数据,.m|

~~~freefem
real Titer = 0.0001; 
real iterMax = 300;

real alphaLow = 2.5*1e-4;
real alphaUp  = 2.5*1e+4;
real qq = 0.1;
real gamma = 0.005;

ofstream object("./Output/Main.obj");

string ModelShapeName = "figure/InitialShape.vtu";
string ModelResultName = "figure/FinialShape.vtu";
string ModelMatLabName = "matlab/Model2_Shape1.m";
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
\nabla u &= [\partial_x u, \partial_y u, \partial_z u] \\
Du &= [\nabla u_x, \nabla u_y, \nabla u_z] \\
\texttt{div} \ u &= \partial_x u_x + \partial_y u_y + \partial_z u_z
\end{align*}
$$

~~~freefem
macro grad(u) [dx(u), dy(u), dz(u)] //
macro Grad(ux, uy, uz) [dx(ux), dy(ux), dz(ux), dx(uy), dy(uy), dz(uy), dx(uz), dy(uz), dz(uz)] //
macro div(ux, uy, uz) ( dx(ux) + dy(uy) + dz(uz) ) //
~~~

## 4. 三维网格和有限元空间
|变量名|作用|
|:---:|:---:|
|Uh|速度场|
|Thetah|下降方向场|
|P1h|Lagrange线性元|
|Phih|$\phi$或$\rho$所在的有限元空间|
|Sth|Uh $\times$ P1h|
|Veh|Thetah $\times$ P1h|
~~~freefem
real d = 1.0;
real r = 0.2;
border Lbottom(t = 0, 1)  {x = t; y = 0; label = 1;}
border Lright(t = 0, 1)   {x = 1; y = t; label = 1;}
border Lup(t = 1, 0)      {x = t; y = 1; label = 1;}
border Lleft(t = 1, 0)    {x = 0; y = t; label = 1;}
border C1(t = 0, 2*pi)    {x = 0.5 + r*cos(t); y = 0.9-r + r*sin(t);}
border C2(t = 0, 2*pi)    {x = 0.5 + r*cos(t); y = 0.5 + r*sin(t);}

int n = 16;
mesh Th2D = buildmesh(Lbottom(n) + Lright(n) + Lup(n) + Lleft(n) + C1(n));
mesh Thbottom2D = buildmesh(Lbottom(n) + Lright(n) + Lup(n) + Lleft(n) + C2(n));
mesh Thup2D = buildmesh(Lbottom(n) + Lright(n) + Lup(n) + Lleft(n));

meshS ThLeft = movemesh23(Th2D, transfo = [x, 0, y]);
meshS Thright = movemesh23(Th2D, transfo = [x, 1, y], orientation = -1);
meshS Thbottom = movemesh23(Thbottom2D, transfo = [y, x, 0]);
meshS Thup = movemesh23(Thup2D, transfo = [y, x, 1], orientation = -1);
meshS Thinfront = movemesh23(Th2D, transfo = [1, x, y]);
meshS Thbehind = movemesh23(Th2D, transfo = [0, x, y], orientation = -1);
meshS ThSurface = ThLeft + Thright + Thbottom + Thup + Thinfront + Thbehind;

checkmesh(ThSurface, rebuildboundary = true);
cout << "Check Over!" << endl;

macro Labelfunc()
(
      3
    - 2 * ( (sqrt((x - 0.5)^2 + (z - 0.9+r)^2) <= r) && ( (abs(y) <= 1e-6) || (abs(y-1) <= 1e-6)  ) )
    - 2 * ( (sqrt((y - 0.5)^2 + (z - 0.9+r)^2) <= r) && ( (abs(x) <= 1e-6) || (abs(x-1) <= 1e-6)  ) )
    - 1 * ( (sqrt((x - 0.5)^2 + (y - 0.5)^2) <= r) && (abs(z) <= 1e-6))
    
)// 
// 1--In, 2--Out, 3--Wall

real dh = 1.0/16.0;
real VolTh = sqrt(2.0) * (dh^3) / 12.0;
mesh3 Th = tetg(ThSurface, switch = "pqaAAQY");
Th = tetgreconstruction(Th, sizeofvolume = VolTh);
Th = change(Th, flabel = Labelfunc);
savemesh(Th, "ThModel2_1.mesh");
savevtk("ThModel2_1.vtu", Th);
savegmsh(Th, "Th");

fespace Uh(Th, [P2, P2, P2]);       
fespace Thetah(Th, [P1, P1, P1]);   
fespace P1h(Th, P1);           
fespace Phih(Th, P1dc3d);       
fespace Sth = Uh * P1h;

P1h h1 = hTriangle;
real hmax = h1[].max; 
real hmin = h1[].min;
object << "maximal size of a mesh triangle h = " << hmax << endl;
object << "minimal size of a mesh triangle h = " << hmin << endl;
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
macro Cy1()
(
    max( sqrt((x-0.5)^2+(z-0.9+r)^2) - r,
         abs(y - 0.5) - 0.5 )
)//
macro Cy2()
(
    max( sqrt((y-0.5)^2+(z-0.9+r)^2) - r,
         abs(x - 0.5) - 0.5 )
)//
macro Cy3()
(
    max( sqrt((x-0.5)^2+(y-0.5)^2) - r,
         abs(z - 0.45) - 0.45 )
)// 
func phi0 = min(Cy1, min(Cy2, Cy3));

P1h phiInit;
distance(Th, phi0, phiInit[]);

Phih phi, Rho;
phi = phiInit;
//Rho = (phi < 0);
Rho = HeaviSide(phi);
plot(Rho, fill = 1, value = 1, cmm = "Rho_0");
savevtk("./figure/phi_0.vtu", Th, phi);
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
Uh [ux, uy, uz], [vx, vy, vz];

macro gvy() 
(
    (abs(y) < 1e-6) * 1.0 + (abs(y-1) < 1e-6) * (-1.0)
)//
macro gvx() 
(
    (abs(x) < 1e-6) * 1.0 + (abs(x-1) < 1e-6) * (-1.0)
)//
func gx = gvx;
func gy = gvy;
func gz = 0;
func fx = 0;
func fy = 0;
func fz = 0;

varf StokesSteady(<[ux, uy, uz], [p]>, <[vx, vy, vz], [q]>) = 
int3d(Th)(
            Grad(ux, uy, uz)' * Grad(vx, vy, vz)
            - p * div(vx, vy, vz)
            - q * div(ux, uy, uz) 
            - p * q * 1e-8
        )
+ on(1, ux = gx, uy = gy, uz = gz)
+ on(3, ux = 0, uy = 0, uz = 0);

varf StokesRHS(<[ux, uy, uz], [p]>, <[vx, vy, vz], [q]>) = 
int3d(Th) (fx * vx + fy * vy + fz * vz)
+ on(1, ux = gx, uy = gy, uz = gz)
+ on(3, ux = 0, uy = 0, uz = 0);

matrix A1m = StokesSteady(Sth, Sth);
real[int] Stokesb = StokesRHS(0, Sth);

macro StokesPDE(ux, uy, uz, Rho)
{
    matrix Am;

    varf StokesRho(<[ux, uy, uz], [p]>, <[vx, vy, vz], [q]>) = 
    int3d(Th)( alpha(Rho) * ([ux, uy, uz]' * [vx, vy, vz]) );
    matrix ARho = StokesRho(Sth, Sth);

    Am = A1m + ARho;

    Mat A(Am);
    cout << A.n << endl;
    /*set(A, sparams = "-pc_type lu");*/

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
Thetah [thetax, thetay, thetaz], [Wx, Wy, Wz];

varf VelocityA([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
int3d(Th)( 
            Grad(thetax, thetay, thetaz)' * Grad(Wx, Wy, Wz)
            + ([thetax, thetay, thetaz]' * [Wx, Wy, Wz])
        )
+ on(1, 2, 3, thetax = 0, thetay = 0, thetaz = 0);

matrix A2m = VelocityA(Thetah, Thetah);

macro VelocityPDE(thetax, thetay, thetaz, ux, uy, uz, Rho)
{
    varf VelocitycRho([thetax, thetay, thetaz], [Wx, Wy, Wz]) = 
    int3d(Th) (div(Wx, Wy, Wz) * Rho)
    + on(1, 2, 3, thetax = 0, thetay = 0, thetaz = 0);

    varf VelocitybRho([thetax, thetay, thetaz], [Wx, Wy, Wz]) =
    int3d(Th)( 
                dalpha(Rho) * (0.5 * ([ux, uy, uz]' * [ux, uy, uz]))  
                * ([Wx, Wy, Wz]' * grad(Rho))
            )
    + on(1, 2, 3, thetax = 0, thetay = 0, thetaz = 0);

    real[int] cRho = VelocitycRho(0, Thetah);
    real[int] bRho = VelocitybRho(0, Thetah);

    matrix LHSm = [[A2m, cRho], [cRho', 0]];
    real[int] RHS = [bRho, 0];
    set(LHSm, solver = CG, eps = 1e-8);
    Mat LHS(LHSm);
    /*set(LHS, sparams = "-pc_type lu");*/

    real[int] sol(RHS.n);
    sol = LHS^-1 * RHS;
    thetax[] = sol(0:sol.n-1);
}//
~~~
计算目标泛函$J(\Omega)$，形状导数$dJ(\Omega;\theta)$以及当前区域体积Vol和$d\texttt{Vol}(\Omega;\theta)$：
~~~freefem
macro GetDrag(ux, uy, uz, Rho)
(
    int3d(Th) ( 0.5 * alpha(Rho) * ([ux, uy, uz]' * [ux, uy, uz]))
    + int3d(Th) ( 0.5 * (Grad(ux, uy, uz)' * Grad(ux, uy, uz)) )
    - int3d(Th) (fx * ux + fy * uy + fz * uz)
) //

macro GetVol(Rho)
(
    int3d(Th) (Rho)
)//

macro GetdJ(ux, uy, uz, thetax, thetay, thetaz, Rho)
(
    int3d(Th)( 
                ( -dalpha(Rho) * (0.5 * ([ux, uy, uz]' * [ux, uy, uz])) ) 
                * [thetax, thetay, thetaz]' * grad(Rho)
            )
)//

macro GetdV(thetax, thetay, thetaz, Rho)
(
    int3d(Th)( 
                [thetax, thetay, thetaz]' * grad(Rho)
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
    N.x * thetax + N.y * thetay + N.z * thetaz
)//

macro LFDG(phiOutput, phiInput, thetax, thetay, thetaz, Tsum)
{
    real t = 0.;
    real dt = Tsum / 5.;

    Phih u, v, phitt, phie;
    P1h OutRho;
    
    varf Adual(u, v) =
    int3d(Th) ( (u / dt + ([thetax, thetay, thetaz]' * grad(u))) * v )
    - intallfaces(Th) (
                        (nElementonB - 1) * v 
                        * (0.5 * abs(nn) - nn / 2) * jump(u) 
                    );
    
    varf DGRHS(u, v) = int3d(Th) ( (u / dt) * v );
         
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
### 迭代时间步长选取
Titer随迭代变化，当前迭代得到$\phi^*$，对于(Titer, $\phi^*$) 

当$J(\phi^*) > J(\phi)[1+\texttt{sign}(J(\phi)) \varepsilon]$时，  
Titer:=Titer/1.5，循环次数＋1

当$J(\phi^*) \leqslant J(\phi)[1+\texttt{sign}(J(\phi)) \varepsilon]$  
Titer:=Titer * 2，迭代次数＋1

~~~freefem
/*
func bool IFchangeTiter(real dragnew, real dragold, int i)
{
    real s = sign(dragold);
    real coe;
    if (i <= 10)
    {
        coe = 0.05;  dTmin = 1e-3;
    }
    else if (10 < i && i <= 120)
    {
        coe = 0.001; dTmin = 1e-6;
    }
    else
    {
        coe = 0.0001; dTmin = 1e-8;
    }
    return ((dragnew > dragold * (1 + s * coe)) && (Titer > dTmin));
}
*/
~~~
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
Uh [uxnew, uynew, uznew];
Phih phitest, phinew, Rhonew;
real[int] drag(iterMax + 1);
real[int] Vol(iterMax + 1);
real dragnew, Volnew;

int iter = 0;
cout << "Yes1" << endl;
StokesPDE(ux, uy, uz, Rho);  
cout << "Yes2" << endl;
VelocityPDE(thetax, thetay, thetaz, ux, uy, uz, Rho);
cout << "Yes3" << endl;
drag[iter] = GetDrag(ux, uy, uz, Rho);
Vol[iter] = GetVol(Rho);
OutPutData(iter, drag, Vol, object);
~~~
### 6.4 迭代求解
~~~freefem
real time = clock();
/*
while ( iter <= iterMax )
{
    cout << "iter = " << iter << endl;
    object << endl;
    object << "iter = " << iter << endl;

    phitest = phi;
    LFDG(phinew, phitest, thetax, thetay, thetaz, Titer);
    if ((iter - 1) % 8 == 0)
        RedistanceConvert(phinew);
    Rhonew = (phinew < 0);
    // Rhonew = HeaviSide(phinew);
    StokesPDE(uxnew, uynew, uznew, Rhonew);
    dragnew = GetDrag(uxnew, uynew, uznew, Rhonew);
    Volnew = GetVol(Rhonew); 

    if ((iter - 1) % 8 == 0)
        cout << "Redistance End!" << endl;
    iter++;
    object << "Titer = " << Titer << endl;
    drag[iter] = dragnew;
    Vol[iter] = Volnew;

    OutPutData(iter, drag, Vol, object);

    phi = phinew;
    [ux, uy, uz] = [uxnew, uynew, uznew];
    Rho = Rhonew;		
    plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter);
    string filename = "./figure/phi_"+iter+".vtu";
    savevtk(filename, Th, phi);
    VelocityPDE(thetax, thetay, thetaz, ux, uy, uz, Rho);

    object << "iter Time = " << clock() - time << endl;
    cout << "iter Time = " << clock() - time << endl;
    time = clock();
    object.flush;
}
*/

while ( iter <= iterMax )
{
    
    cout << "iter = " << iter << endl;
    object << endl;
    object << "iter = " << iter << endl;

    phitest = phi;
    LFDG(phi, phitest, thetax, thetay, thetaz, Titer);
    if ((iter - 1) % 10 == 0)
    {
        RedistanceConvert(phi);
        cout << "Redistance End!" << endl;
    }
    //Rho = (phi < 0);
    Rho = HeaviSide(phi);
    StokesPDE(ux, uy, uz, Rho);
    iter++;
    object << "Titer = " << Titer << endl;
    drag[iter] = GetDrag(ux, uy, uz, Rho);
    Vol[iter] = GetVol(Rho); 
   
    OutPutData(iter, drag, Vol, object);
	
    plot(Rho, fill = 1, value = 1, cmm = "Rho_"+iter, wait = false);
    string filename = "./figure/phi_"+iter+".vtu";
    if (iter % 10 == 0)
    {
        savevtk(filename, Th, phi, [ux, uy, uz]);
    }
    VelocityPDE(thetax, thetay, thetaz, ux, uy, uz, Rho);

    object << "iter Time = " << clock() - time << endl;
    cout << "iter Time = " << clock() - time << endl;
    time = clock();
    object.flush;
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