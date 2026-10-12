param(
    [string]$FreeFem = 'D:\MathTool\FreeFem++\FreeFem++.exe',
    [int[]]$MeshScales = @(4, 5, 6, 8),
    [double]$BaseT = 0.075,
    [double]$BaseGamma = 0.005,
    [int]$BaseMeshScale = 4,
    [int]$Iterations = 300,
    [switch]$Parallel,
    [int]$MaxParallel = 3,
    [int]$SingleScale = 0,
    [string]$Tag = '',
    [switch]$FixedT,
    [string]$Scales = '',
    [string]$TList = '',
    [double]$LevelT = -1
)

# `powershell -File script.ps1 -MeshScales 3,4,5,6` cannot deliver an int[]: the argument
# binder reads that comma list as the single integer 3456, because commas are treated as
# digit-group separators when converting a string to an integer. Accept an explicit string
# list as well, so -File invocations (the only form the execution policy permits here) work.
if (-not [string]::IsNullOrWhiteSpace($Scales)) {
    $MeshScales = @($Scales -split '[,\s;]+' |
        Where-Object { $_ -ne '' } |
        ForEach-Object { [int]$_ })
}

# Optional explicit per-level T, one value per mesh scale. This is what allows a true
# T proportional to h design: T_j = c * h_max,j using the measured h_max of each mesh.
$explicitT = @()
if (-not [string]::IsNullOrWhiteSpace($TList)) {
    $explicitT = @($TList -split '[,\s;]+' |
        Where-Object { $_ -ne '' } |
        ForEach-Object { [double]$_ })
    if ($explicitT.Count -ne $MeshScales.Count) {
        throw ("TList must hold one value per mesh scale: got {0} T values for {1} scales." -f `
            $explicitT.Count, $MeshScales.Count)
    }
}

$ErrorActionPreference = 'Stop'
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$mainFile = Join-Path $projectRoot 'Main.edp'
$outputDirectory = Join-Path $projectRoot 'Output'
# Tag isolates one study's artifacts AND its temp root, so that a new experiment can
# never overwrite or delete a study that is still running. An empty Tag reproduces the
# original coupled-mesh (T = O(h)) study exactly.
$tagPrefix = if ([string]::IsNullOrWhiteSpace($Tag)) { '' } else { $Tag }
$artifactPrefix = if ($tagPrefix) { 'CoupledMesh' + $tagPrefix } else { 'CoupledMesh' }
$temporaryRoot = Join-Path $projectRoot $(if ($tagPrefix) {
        '.coupled-mesh-' + $tagPrefix + '-work'
    }
    else {
        '.coupled-mesh-study-work'
    })
$finalCsv = Join-Path $outputDirectory ($artifactPrefix + 'Study.csv')

if (-not (Test-Path -LiteralPath $FreeFem -PathType Leaf)) {
    throw "FreeFem++ executable not found: $FreeFem"
}
if (-not (Test-Path -LiteralPath $mainFile -PathType Leaf)) {
    throw "Model source not found: $mainFile"
}
if ($BaseT -le 0 -or $BaseGamma -le 0 -or $BaseMeshScale -lt 1) {
    throw 'BaseT, BaseGamma, and BaseMeshScale must be positive.'
}
if ($Iterations -lt 1) {
    throw 'Iterations must be a positive integer.'
}
if ($MaxParallel -lt 1) {
    throw 'MaxParallel must be positive.'
}

# A per-level CSV is a usable checkpoint only when it holds exactly $Iterations rows.
function Get-CompletedLevelRows {
    param(
        [string]$CsvPath,
        [int]$Expected
    )
    if (-not (Test-Path -LiteralPath $CsvPath -PathType Leaf)) {
        return $null
    }
    $rows = @(Import-Csv -LiteralPath $CsvPath)
    if ($rows.Count -ne $Expected) {
        Write-Output ("Ignoring incomplete checkpoint {0}: {1} of {2} rows" -f `
            (Split-Path -Leaf $CsvPath), $rows.Count, $Expected)
        return $null
    }
    return ,$rows
}

function Get-LevelT {
    param([int]$Scale)
    if ($LevelT -gt 0) {
        # Worker mode: the driver already resolved this level's T and passed it in.
        return $LevelT
    }
    if ($explicitT.Count -gt 0) {
        $idx = [array]::IndexOf($MeshScales, $Scale)
        if ($idx -lt 0) {
            throw "Mesh scale $Scale is not present in the supplied scale list."
        }
        return $explicitT[$idx]
    }
    if ($FixedT) {
        # Fixed-T / fixed-gamma study: only h refines, T is held constant.
        return $BaseT
    }
    return $BaseT * ([double]$BaseMeshScale / [double]$Scale)
}

function Get-PerLevelCsv {
    param([int]$Scale)
    return (Join-Path $outputDirectory ("{0}_n{1}.csv" -f $artifactPrefix, $Scale))
}

# Run exactly one mesh level and store its checkpoint CSV. Used both serially and
# as the parallel worker entry point.
function Invoke-MeshLevel {
    param([int]$Scale)

    if ($Scale -lt 1) {
        throw "Invalid mesh scale: $Scale"
    }
    New-Item -ItemType Directory -Force -Path $outputDirectory | Out-Null

    $caseDirectory = Join-Path $temporaryRoot ("mesh-{0}" -f $Scale)
    $caseOutput = Join-Path $caseDirectory 'Output'
    $caseMatlab = Join-Path $caseDirectory 'matlab'
    New-Item -ItemType Directory -Force -Path $caseOutput, $caseMatlab | Out-Null

    $timeData = Join-Path $caseOutput 'StepDiagnostics.csv'
    $orderData = Join-Path $caseOutput 'UnusedStepOrder.csv'
    $objectData = Join-Path $caseOutput 'Main.obj'
    $matlabData = Join-Path $caseMatlab 'Model.m'
    $caseT = Get-LevelT -Scale $Scale

    Push-Location -LiteralPath $caseDirectory
    try {
        & $FreeFem -v 0 $mainFile `
            -meshScale $Scale `
            -gamma $BaseGamma.ToString('R', $culture) `
            -Titer $caseT.ToString('R', $culture) `
            -iterMax $Iterations `
            -writeVisual 0 `
            -silent 1 `
            -enableShadow 0 `
            -object $objectData `
            -matlab $matlabData `
            -timeData $timeData `
            -stepOrderData $orderData `
            -options_left no
        if ($LASTEXITCODE -ne 0) {
            throw "FreeFem++ failed for mesh scale $Scale with exit code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }

    $rows = @(Import-Csv -LiteralPath $timeData)
    if ($rows.Count -ne $Iterations) {
        throw "Expected $Iterations diagnostic rows for mesh scale $Scale; found $($rows.Count)."
    }
    $perLevelCsv = Get-PerLevelCsv -Scale $Scale
    $rows | Export-Csv -LiteralPath $perLevelCsv -NoTypeInformation -Encoding utf8
    Write-Output "Completed mesh scale ${Scale}: $($rows.Count) rows -> $perLevelCsv"
}

# ---------------------------------------------------------------- worker mode
# Spawned by the parallel driver. Runs one level inside the shared temp root
# that the driver created, and never touches the merged outputs.
if ($SingleScale -gt 0) {
    $existing = Get-CompletedLevelRows -CsvPath (Get-PerLevelCsv -Scale $SingleScale) `
        -Expected $Iterations
    if ($null -ne $existing) {
        Write-Output "Mesh scale $SingleScale is already complete; worker has nothing to do."
        exit 0
    }
    Invoke-MeshLevel -Scale $SingleScale
    exit 0
}

# ---------------------------------------------------------------- driver mode
if ($MeshScales.Count -lt 1) {
    throw 'At least one mesh scale is required.'
}
# Guard against the -File argument trap: 'powershell -File script.ps1 -MeshScales 3,4,5,6'
# delivers the comma list as the single integer 3456 instead of an array.
if ($MeshScales.Count -eq 1 -and [int]$MeshScales[0] -gt 1000) {
    throw ('Mesh scales arrived as the single value {0}, which means a comma list was ' +
        'collapsed by the -File argument binder. Invoke with the call operator instead: ' +
        '& .\RunCoupledMeshStudy.ps1 -MeshScales 3,4,5,6 ...') -f $MeshScales[0]
}
if (($MeshScales | Sort-Object -Unique).Count -ne $MeshScales.Count) {
    throw 'Mesh scales must be distinct.'
}

New-Item -ItemType Directory -Force -Path $outputDirectory | Out-Null
if (Test-Path -LiteralPath $temporaryRoot) {
    $resolved = (Resolve-Path -LiteralPath $temporaryRoot).Path
    if (-not $resolved.StartsWith($projectRoot + [IO.Path]::DirectorySeparatorChar,
            [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove unexpected path: $resolved"
    }
    Remove-Item -LiteralPath $resolved -Recurse -Force
}
New-Item -ItemType Directory -Path $temporaryRoot | Out-Null

try {
    $pending = @()
    foreach ($scale in $MeshScales) {
        $perLevelCsv = Get-PerLevelCsv -Scale $scale
        $rows = Get-CompletedLevelRows -CsvPath $perLevelCsv -Expected $Iterations
        if ($null -ne $rows) {
            Write-Output ("Reusing completed mesh scale {0}: {1} rows <- {2}" -f `
                $scale, $rows.Count, $perLevelCsv)
        }
        else {
            $pending += $scale
        }
    }

    if ($pending.Count -gt 0) {
        $useParallel = $Parallel -and $pending.Count -gt 1
        Write-Output ("Pending mesh scales: {0} (mode: {1})" -f `
            ($pending -join ','), $(if ($useParallel) { 'parallel' } else { 'serial' }))

        if ($useParallel) {
            $scriptPath = $MyInvocation.MyCommand.Path
            $queue = New-Object System.Collections.Queue
            foreach ($scale in $pending) { $queue.Enqueue($scale) }
            $running = @()
            $failures = @()

            while ($queue.Count -gt 0 -or $running.Count -gt 0) {
                while ($running.Count -lt $MaxParallel -and $queue.Count -gt 0) {
                    $scale = [int]$queue.Dequeue()
                    $logFile = Join-Path $temporaryRoot ("mesh-{0}.log" -f $scale)
                    $errFile = Join-Path $temporaryRoot ("mesh-{0}.err.log" -f $scale)
                    $childArgs = @(
                        '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', $scriptPath,
                        '-SingleScale', $scale,
                        '-Iterations', $Iterations,
                        '-BaseMeshScale', $BaseMeshScale,
                        '-BaseGamma', $BaseGamma.ToString('R', $culture),
                        '-BaseT', $BaseT.ToString('R', $culture),
                        '-MeshScales', ($MeshScales -join ','),
                        '-FreeFem', $FreeFem
                    )
                    if (-not [string]::IsNullOrWhiteSpace($Tag)) {
                        $childArgs += @('-Tag', $Tag)
                    }
                    if ($FixedT) {
                        $childArgs += '-FixedT'
                    }
                    # A worker receives no valid MeshScales (the -File binder collapses the
                    # list), so it cannot resolve T itself: the driver resolves and passes it.
                    $childArgs += @('-LevelT', (Get-LevelT -Scale $scale).ToString('R', $culture))
                    Write-Output "Launching mesh scale $scale (log: $logFile)"
                    $proc = Start-Process -FilePath 'powershell' -ArgumentList $childArgs `
                        -NoNewWindow -PassThru `
                        -RedirectStandardOutput $logFile -RedirectStandardError $errFile
                    $running += [pscustomobject]@{
                        scale = $scale; proc = $proc; log = $logFile; err = $errFile
                    }
                }

                Start-Sleep -Seconds 5
                $stillRunning = @()
                foreach ($w in $running) {
                    if ($w.proc.HasExited) {
                        $code = $w.proc.ExitCode
                        if ($null -ne $code -and $code -ne 0) {
                            $tail = if (Test-Path $w.err) {
                                (Get-Content -LiteralPath $w.err -Tail 20) -join "`n"
                            } else { '(no stderr)' }
                            $failures += "mesh scale $($w.scale) exited $code`n$tail"
                        }
                        else {
                            # Start-Process can report a null ExitCode even on success, so the
                            # per-level artifact check in the merge phase is the real gate.
                            Write-Output "Worker finished mesh scale $($w.scale)"
                        }
                    }
                    else {
                        $stillRunning += $w
                    }
                }
                $running = $stillRunning
            }

            if ($failures.Count -gt 0) {
                throw ($failures -join "`n`n")
            }
        }
        else {
            foreach ($scale in $pending) {
                Invoke-MeshLevel -Scale $scale
            }
        }
    }

    # ------------------------------------------------------------- merge phase
    $allRows = @()
    $levelMetadata = @()
    $level = 0
    foreach ($scale in $MeshScales) {
        $perLevelCsv = Get-PerLevelCsv -Scale $scale
        $rows = Get-CompletedLevelRows -CsvPath $perLevelCsv -Expected $Iterations
        if ($null -eq $rows) {
            throw "Mesh scale $scale is still incomplete after the study run."
        }
        foreach ($row in $rows) {
            $row | Add-Member -NotePropertyName coupled_level -NotePropertyValue $level -Force
            $row | Add-Member -NotePropertyName transport_relative_drift `
                -NotePropertyValue ([double]$row.transport_abs_drift / [double]$row.volume_before) -Force
            $allRows += $row
        }
        # Re-export so every per-level checkpoint carries identical columns.
        $rows | Export-Csv -LiteralPath $perLevelCsv -NoTypeInformation -Encoding utf8

        $first = $rows[0]
        $levelMetadata += [pscustomobject]@{
            coupled_level = $level
            mesh_scale    = [int]$first.mesh_scale
            hmax          = [double]$first.hmax
            hmin          = [double]$first.hmin
            gamma         = [double]$first.gamma
            T             = [double]$first.T
            T_over_hmax   = [double]$first.T / [double]$first.hmax
            rows          = $rows.Count
        }
        $level++
    }

    $temporaryCsv = Join-Path $temporaryRoot 'CoupledMeshStudy.csv'
    $allRows | Export-Csv -LiteralPath $temporaryCsv -NoTypeInformation -Encoding utf8
    Copy-Item -LiteralPath $temporaryCsv -Destination $finalCsv -Force
    $metadataCsv = Join-Path $outputDirectory ($artifactPrefix + 'Levels.csv')
    $levelMetadata | Export-Csv -LiteralPath $metadataCsv -NoTypeInformation -Encoding utf8

    $orderedMetadata = @($levelMetadata | Sort-Object hmax -Descending)
    $ratioH = $orderedMetadata[-1].hmax / $orderedMetadata[0].hmax
    Write-Output "Saved $($allRows.Count) coupled mesh diagnostic rows: $finalCsv"
    Write-Output "Saved level metadata: $metadataCsv"
    Write-Output ("Actual finest/coarsest hmax ratio: {0:R}" -f $ratioH)
}
finally {
    if (Test-Path -LiteralPath $temporaryRoot) {
        $resolved = (Resolve-Path -LiteralPath $temporaryRoot).Path
        if (-not $resolved.StartsWith($projectRoot + [IO.Path]::DirectorySeparatorChar,
                [StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove unexpected path: $resolved"
        }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    }
}
