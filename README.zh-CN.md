# Volume-Controlled Level-Set Method / 体积控制水平集方法

[English README](README.md)

本项目收集基于 FreeFEM 的体积控制水平集形状与拓扑优化算例，包含二维和三维线弹性、弹性特征值问题及 Stokes 流动，并提供 MATLAB 后处理、Gmsh 网格定义以及时间步和体积误差诊断。

This repository contains FreeFEM research examples for volume-controlled
level-set shape and topology optimization, with MATLAB post-processing,
Gmsh geometry definitions, and transport/reinitialization diagnostics.

## Directory Layout / 目录结构

| Directory / 目录 | Contents / 内容 |
| --- | --- |
| `Main Code/2DElastic/` | 2D linear elasticity / 二维线弹性 |
| `Main Code/3DElastic/` | 3D linear elasticity / 三维线弹性 |
| `Main Code/2DElastic eigenvalue/` | 2D elastic eigenvalues / 二维弹性特征值 |
| `Main Code/3DElastic eigenvalue/` | 3D elastic eigenvalues / 三维弹性特征值 |
| `Main Code/2DStokes/` | 2D Stokes flow / 二维 Stokes 流动 |
| `Main Code/3DStokes/` | 3D Stokes flow / 三维 Stokes 流动 |
| `Estimate For Time/` | Time/mesh refinement and volume diagnostics / 时间步、网格加密和体积诊断 |
| `GenMesh_3D_Model1/` | First 3D geometry family / 第一类三维几何及网格 |
| `GenMesh_3D_Model2/` | Second 3D geometry family / 第二类三维几何及网格 |

保留原有算例路径和中文说明。各算例通常包含 `Main.md`、`matlab/`、图像及输出数据。

The original example paths and notes are retained. Each example typically
contains `Main.md`, MATLAB scripts or histories, figures, and output data.

## Requirements / 环境依赖

- **FreeFEM 4.15+**: supports direct execution of Markdown examples / 支持直接运行 Markdown 算例。
- **MATLAB**: history analysis and plotting / 历史数据分析和绘图。
- **Gmsh**: mesh generation / 网格生成。
- **ParaView**: optional VTK/VTU visualization / 可选的结果可视化工具。
- **PowerShell**: optional refinement-study drivers / 可选的加密研究驱动脚本环境。

FreeFEM 插件因算例而异，常见依赖包括 `distance`、`iovtk`、`gmsh`、`TetGen`/`tetgen`、`PETSc`，部分算例还使用 `msh3` 或 `MUMPS_seq`。MATLAB 脚本中的 `exportgraphics`、`prctile` 等函数需要对应版本或工具箱支持。

Plugin and solver availability depends on the installation. See the
[FreeFEM changelog](https://github.com/FreeFem/FreeFem-sources/blob/develop/CHANGELOG.md)
for Markdown execution support introduced in version 4.15.

## Execution Reference / 运行方式参考

以下命令用于说明作者的运行流程。许可证仅允许下载和查看；其他使用（包括运行代码）需要事先取得书面许可。

The commands below document the authors' workflow. The license grants download
and inspection rights only; execution and other uses require prior written permission.

在目标算例目录内运行，确保相对网格路径正确：

Run from the selected example directory so that relative paths resolve:

```powershell
New-Item -ItemType Directory -Force -Path Output, figure, matlab | Out-Null
FreeFem++ Main.md
```

无需手动从 Markdown 提取 `.edp`。已有独立脚本的算例也可使用其 `.edp` 文件。`Main.md` 和 `.edp` 不保证完全一致，应核对实际运行文件。完整优化可能耗时较长。

No manual Markdown extraction is required. Standalone `.edp` scripts are also
available in the time study and selected eigenvalue examples. Markdown and
standalone files are not guaranteed to be identical.

## Time and Volume Study / 时间步与体积诊断

`Estimate For Time/Main.edp` 提供 `-Titer`、`-iterMax`、`-gamma`、`-meshScale`、`-enableShadow` 等参数。`RunFullStepStudy.ps1` 和 `RunCoupledMeshStudy.ps1` 用于批量研究；其默认 FreeFEM 路径属于作者本机，其他环境需通过 `-FreeFem` 指定。

Use distinct output paths for independent runs. The study drivers support a
custom FreeFEM executable path through `-FreeFem`.

CSV 数据分别记录 `volume_before`、`volume_after_transport` 和 `volume_after_step`；`reinit_abs_jump` 直接表示重初始化引起的绝对体积跳变，不能通过两个绝对误差相减代替。

The diagnostics distinguish transport changes from reinitialization jumps.
The numerical outputs do not guarantee exact discrete volume conservation.

## Meshes and Post-Processing / 网格与后处理

三维算例需保留源码要求的 `.msh` 文件名和物理边界标签。重新生成网格前应检查 `.geo` 中的参数和 `Save` 语句，不能仅凭文件名判断实际网格尺寸。

Preserve mesh filenames and physical labels expected by each example.
Review geometry parameters and output statements before regeneration.

MATLAB 绘图脚本可能依赖 `Num`、`Volume`、`compliance` 等工作区变量，应先加载对应历史数据。部分 `.m` 文件是导出的数值历史。绘图可能覆盖同名 PDF，保留参考结果时应使用独立副本。

Load the corresponding history before plotting. Some `.m` files contain data
rather than reusable functions. Preserve reference outputs before rerunning.

## Data Policy / 数据管理

`.gitignore` 排除 `.vtu`、`.vtk` 可视化结果、`.obj` 文本日志、临时研究目录和本地编辑器配置；保留源码说明、MATLAB 文件、CSV、PNG/PDF 和普通输入网格。

The oversized `GenMesh_3D_Model1/Msh/Th003_002.msh` is excluded; its geometry
definition remains available. Full visualization sequences and this mesh can
be distributed separately or managed with Git LFS. Local files are not deleted.

当前参考结果尚未通过干净环境完整复现验证；复现需要匹配求解器环境及原始迭代设置。

A clean checkout has not been verified to reproduce all supplied figures.

## License and Citation / 许可证与引用

本项目原创材料采用自定义 [Download and Inspection Only License](LICENSE)：允许下载和查看未修改副本；禁止修改、改编、衍生作品和再分发，其他使用需要事先获得书面许可。GitHub 平台条款和法律要求的权利除外。第三方材料遵循各自许可证。

This is a restrictive source-available license, not an open-source license.
The copyright holders may continue maintaining their own materials.
Publication citation metadata has not yet been specified / 论文引用信息待补充。
