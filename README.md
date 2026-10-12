# Volume-Controlled Level-Set Method

[Chinese-English README](README.zh-CN.md)

FreeFEM research examples for volume-controlled level-set shape and topology
optimization, with MATLAB post-processing and Gmsh geometry definitions.
The examples cover two- and three-dimensional linear elasticity, elastic
eigenvalue problems, and Stokes flow. A separate study directory contains
diagnostics for volume changes during level-set transport and reinitialization.

## Repository Layout

| Directory | Contents |
| --- | --- |
| `Main Code/2DElastic/` | Two-dimensional linear elasticity examples and variants |
| `Main Code/3DElastic/` | Three-dimensional linear elasticity examples and shared MATLAB post-processing |
| `Main Code/2DElastic eigenvalue/` | Two-dimensional elastic eigenvalue optimization examples |
| `Main Code/3DElastic eigenvalue/` | Three-dimensional elastic eigenvalue optimization example |
| `Main Code/2DStokes/` | Two-dimensional Stokes flow examples |
| `Main Code/3DStokes/` | Three-dimensional Stokes flow examples |
| `Estimate For Time/` | Transport, reinitialization, and time/mesh refinement diagnostics |
| `GenMesh_3D_Model1/` | Gmsh geometry definitions and meshes for the first 3D geometry family |
| `GenMesh_3D_Model2/` | Gmsh geometry definition and mesh for the second 3D geometry family |

Individual example directories typically contain `Main.md`, a `matlab/`
directory, figures, and output data. Some directory names and explanatory notes
in the existing examples are in Chinese; this README provides an English entry
point while retaining the original paths.

## Numerical Workflow

The examples combine finite element state solves, computation of a design
velocity, discontinuous Galerkin level-set transport, and distance-based
reinitialization. Objective and volume histories are exported for analysis.
Boundary conditions, initial shapes, discretizations, solver choices, and
optimization parameters are specified separately in each example.

Volume control is assessed through numerical diagnostics. The supplied outputs
should not be interpreted as a guarantee of exact discrete volume conservation
for every transport or reinitialization step.

## Requirements

- **FreeFEM 4.15 or later** for direct execution of the Markdown examples,
  with the plugins required by the selected example. Common loads
  include `distance`, `iovtk`, `gmsh`, `TetGen`/`tetgen`, and `PETSc`. Some examples
  also load `msh3` or `MUMPS_seq`. Eigenvalue examples use `EigenValue`.
- **MATLAB** for the supplied `.m` scripts. Some scripts use functions such as
  `exportgraphics` and `prctile`; check the selected script against your MATLAB
  release and installed toolboxes.
- **Gmsh** to regenerate or modify the supplied 3D meshes.
- **ParaView**, optionally, to inspect VTK/VTU outputs and saved `.pvsm` states.
- **PowerShell**, optionally, to use the refinement-study drivers.

The Markdown execution feature was introduced in FreeFEM 4.15; see the
[official changelog](https://github.com/FreeFem/FreeFem-sources/blob/develop/CHANGELOG.md).
Minimum versions for all other features and dependencies have not been
established. Plugin availability and solver support depend on the installation.

## Running an Example

The commands below document the authors' workflow. The license grants download
and inspection rights only; execution and other uses require prior written permission.

Run FreeFEM from the example directory so that relative mesh and output paths
resolve correctly. Create `Output/`, `figure/`, and `matlab/` if they are absent.

Most examples store their FreeFEM source in fenced `freefem` code blocks inside
`Main.md`. FreeFEM 4.15 and later can execute Markdown scripts directly,
running the FreeFEM code blocks while ignoring the surrounding documentation.
From the selected example directory, run:

```powershell
New-Item -ItemType Directory -Force -Path Output, figure, matlab | Out-Null
FreeFem++ Main.md
```

No manual extraction or conversion to `.edp` is required with a compatible
FreeFEM installation. Review the example's parameters before execution; a full
optimization run can be computationally expensive.

Standalone FreeFEM scripts are currently available in:

- `Estimate For Time/Main.edp`
- The Model1/Shape1 example under `Main Code/2DElastic eigenvalue/` (`main.edp`)
- The Model1/Shape1 example under `Main Code/3DElastic eigenvalue/` (`Main.edp`)

These are separate maintained files; do not assume that every `Main.md` and
`.edp` pair is identical. Review the actual script before execution.

### Time and Volume Diagnostics

The following PowerShell commands run a short diagnostic example. They assume
that `FreeFem++` is available on `PATH` and are issued from the repository root:

```powershell
Set-Location 'Estimate For Time'
New-Item -ItemType Directory -Force -Path Output, figure, matlab | Out-Null
FreeFem++ Main.edp -iterMax 5 -writeVisual 0 -silent 1 `
  -object './Output/QuickCheck.obj' `
  -matlab './matlab/QuickCheck.m' `
  -timeData './Output/QuickCheckDiagnostics.csv' `
  -stepOrderData './Output/QuickCheckStepOrder.csv'
```

This is an example invocation, not a recorded validation run. Five iterations
are intended only as an initial execution check, not as a convergence study.

Parameters exposed by this script include `-Titer`, `-iterMax`, `-gamma`,
`-meshScale`, `-enableShadow`, `-shadowStride`, and `-shadowLevels`. Shadow probes
are disabled by default. Use distinct output paths for independent runs.

`RunFullStepStudy.ps1` and `RunCoupledMeshStudy.ps1` drive longer studies.
Their default FreeFEM executable path is specific to the original Windows
machine; supply `-FreeFem` with your own executable path. Review their parameters
and output filenames before starting a study.

Diagnostic CSV fields distinguish `volume_before`, `volume_after_transport`,
and `volume_after_step`. In particular, `reinit_abs_jump` measures the absolute
volume jump due to reinitialization. It is not obtained by subtracting two
absolute volume errors. MATLAB scripts in `Estimate For Time/matlab/` produce
summary tables and figures from these diagnostics.

### Three-Dimensional Meshes

Several 3D examples read local `.msh` files with `gmshload3`. Preserve the mesh
filenames and physical boundary labels expected by the selected source.
The supplied `.geo` files define geometry and mesh-generation settings.
Review their mesh-size parameters and `Save` statements before regeneration;
filenames alone should not be used to infer the actual mesh resolution.

### MATLAB Post-Processing

Run the selected example's history script in MATLAB before its plotting script
when the plot expects workspace variables such as `Num`, `Volume`, or
`compliance`. Some `.m` files contain exported numerical histories rather than
general-purpose functions. Plot scripts may write PDF figures to the current
working directory; use a separate copy when preserving reference outputs.

## Data and Output Policy

The default `.gitignore` excludes generated `.vtu` and `.vtk` files, text run
logs stored as `.obj`, temporary study directories, and local editor settings.
It retains source documents, MATLAB scripts and histories, CSV diagnostics,
PNG/PDF figures, and ordinary input meshes.

The oversized `GenMesh_3D_Model1/Msh/Th003_002.msh` is also excluded. Its geometry
definition remains available. Distribute this mesh and full visualization
sequences separately, or configure Git LFS before adding them to Git history.
Excluded files remain on the local machine.

Existing output data are reference artifacts. A clean checkout has not been
verified to reproduce all supplied figures, and complete reproduction may
require matching solver environments and running the original iteration counts.

## Citation and License

The original project materials are provided under the custom
[Download and Inspection Only License](LICENSE). Downloading and inspecting
unmodified copies is permitted. Modification, derivative works, redistribution,
and other uses require prior written permission, except for rights required by
applicable law or GitHub's Terms of Service. This is not an open-source license.
Separately licensed third-party materials retain their own terms.

Publication citation metadata has not yet been specified.
