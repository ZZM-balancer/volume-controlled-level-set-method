# Volume-Controlled Level-Set Method

[Chinese-English README](README.zh-CN.md)

FreeFEM examples for volume-controlled level-set shape and topology
optimization in two and three dimensions, covering linear elasticity,
elastic eigenvalue problems, and Stokes flow.

## Contents

| Directory | Description |
| --- | --- |
| `Main Code/2DElastic/` | 2D linear elasticity |
| `Main Code/3DElastic/` | 3D linear elasticity |
| `Main Code/2DElastic eigenvalue/` | 2D elastic eigenvalue optimization |
| `Main Code/3DElastic eigenvalue/` | 3D elastic eigenvalue optimization |
| `Main Code/2DStokes/` | 2D Stokes flow |
| `Main Code/3DStokes/` | 3D Stokes flow |
| `Estimate For Time/` | Time/mesh refinement and volume-error diagnostics |
| `GenMesh_3D_Model1/`, `GenMesh_3D_Model2/` | Gmsh geometry definitions and meshes |

Each example includes documented FreeFEM code in `Main.md`, together with
available MATLAB histories, plotting scripts, and reference figures.

## Requirements

- **FreeFEM 4.15+** for direct execution of Markdown scripts. Required plugins
  vary by example and include `distance`, `iovtk`, `gmsh`, `PETSc`, and others.
- **MATLAB** for post-processing.
- **Gmsh** for mesh generation and **ParaView** for visualization, when needed.

## Running an Example

Open a terminal in the selected example directory and run:

```powershell
New-Item -ItemType Directory -Force -Path Output, figure, matlab | Out-Null
FreeFem++ Main.md
```

Review parameters and mesh paths in `Main.md` before running. FreeFEM executes
the embedded code blocks directly. Examples with standalone `.edp` files can
also be run using those scripts.

For MATLAB plots, load the corresponding history file before running the
plotting script. Time and volume studies are located in `Estimate For Time/`;
the PowerShell drivers accept a custom executable path through `-FreeFem`.

## Output Files

The repository retains scripts, CSV diagnostics, reference figures, and ordinary
input meshes. Bulk VTK/VTU outputs and the oversized `Th003_002.msh` mesh are
excluded from Git; the corresponding Gmsh geometry definition is retained.
