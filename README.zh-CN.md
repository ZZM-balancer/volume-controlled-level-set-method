# Volume-Controlled Level-Set Method / 体积控制水平集方法

[English README](README.md)

基于 FreeFEM 的二维和三维体积控制水平集形状与拓扑优化算例，涵盖线弹性、弹性特征值问题和 Stokes 流动。

FreeFEM examples for volume-controlled level-set shape and topology
optimization in linear elasticity, elastic eigenvalue problems, and Stokes flow.

## Contents / 目录结构

| Directory / 目录 | Description / 内容 |
| --- | --- |
| `Main Code/2DElastic/` | 2D linear elasticity / 二维线弹性 |
| `Main Code/3DElastic/` | 3D linear elasticity / 三维线弹性 |
| `Main Code/2DElastic eigenvalue/` | 2D elastic eigenvalues / 二维弹性特征值 |
| `Main Code/3DElastic eigenvalue/` | 3D elastic eigenvalues / 三维弹性特征值 |
| `Main Code/2DStokes/` | 2D Stokes flow / 二维 Stokes 流动 |
| `Main Code/3DStokes/` | 3D Stokes flow / 三维 Stokes 流动 |
| `Estimate For Time/` | Time/mesh refinement and volume diagnostics / 时间步、网格加密和体积诊断 |
| `GenMesh_3D_Model1/`, `GenMesh_3D_Model2/` | Gmsh geometries and meshes / 几何定义及网格 |

各算例包含 `Main.md` 中的 FreeFEM 代码及说明，并附有已有的 MATLAB 历史数据、绘图脚本和参考图像。

## Requirements / 环境依赖

- **FreeFEM 4.15+**: direct Markdown execution / 直接运行 Markdown 算例。
- **MATLAB**: post-processing / 数据分析及绘图。
- **Gmsh** and **ParaView**: mesh generation and visualization when needed / 按需生成网格和查看结果。

FreeFEM 插件因算例而异，常见依赖包括 `distance`、`iovtk`、`gmsh` 和 `PETSc` 等，应以各算例的 `load` 语句为准。

## Running an Example / 运行算例

在目标算例目录打开终端并执行 / Run from the selected example directory:

```powershell
New-Item -ItemType Directory -Force -Path Output, figure, matlab | Out-Null
FreeFem++ Main.md
```

运行前检查 `Main.md` 中的参数和网格路径。FreeFEM 会直接执行其中的代码块；已有独立 `.edp` 文件的算例也可使用相应脚本。

Review parameters and mesh paths before execution. Markdown code blocks run
directly; standalone `.edp` scripts are also available for selected examples.

MATLAB 绘图前先加载对应历史数据文件。时间步与体积诊断位于 `Estimate For Time/`，其中 PowerShell 驱动可通过 `-FreeFem` 指定本机可执行文件路径。

## Output Files / 输出文件

保留脚本、CSV 诊断、参考图像及普通输入网格。批量 VTK/VTU 结果和超大网格 `Th003_002.msh` 不纳入 Git，对应 Gmsh 几何定义仍保留。

Scripts, diagnostics, reference figures, and ordinary input meshes are retained;
bulk visualization outputs and the oversized mesh are excluded from Git.
