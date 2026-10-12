param(
    [string]$FreeFem = 'D:\MathTool\FreeFem++\FreeFem++.exe',
    [int]$Iterations = 300,
    [double]$T = 0.075,
    [int]$ShadowStride = 10,
    [int]$ShadowLevels = 4
)

$ErrorActionPreference = 'Stop'
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$mainFile = Join-Path $projectRoot 'Main.edp'
$outputDirectory = Join-Path $projectRoot 'Output'
$temporaryRoot = Join-Path $projectRoot '.full-step-study-work'
$stepData = Join-Path $outputDirectory 'StepDiagnostics.csv'
$orderData = Join-Path $outputDirectory 'StepOrder.csv'

if (-not (Test-Path -LiteralPath $FreeFem -PathType Leaf)) {
    throw "FreeFem++ executable not found: $FreeFem"
}
if (-not (Test-Path -LiteralPath $mainFile -PathType Leaf)) {
    throw "Model source not found: $mainFile"
}
if ($Iterations -lt 1) {
    throw 'Iterations must be positive.'
}
if ($T -le 0.0) {
    throw 'T must be positive.'
}
if ($ShadowStride -lt 1) {
    throw 'ShadowStride must be positive.'
}
if ($ShadowLevels -lt 3) {
    throw 'At least three shadow time levels are required.'
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
New-Item -ItemType Directory -Path (Join-Path $temporaryRoot 'Output') | Out-Null
New-Item -ItemType Directory -Path (Join-Path $temporaryRoot 'matlab') | Out-Null
New-Item -ItemType Directory -Path (Join-Path $temporaryRoot 'figure') | Out-Null

try {
    Push-Location -LiteralPath $temporaryRoot
    try {
        & $FreeFem -v 0 $mainFile `
            -Titer $T.ToString('R', $culture) `
            -iterMax $Iterations `
            -writeVisual 0 `
            -silent 1 `
            -enableShadow 1 `
            -shadowStride $ShadowStride `
            -shadowLevels $ShadowLevels `
            -object '.\Output\Main.obj' `
            -matlab '.\matlab\Model1_Shape1.m' `
            -timeData $stepData `
            -stepOrderData $orderData `
            -options_left no
        if ($LASTEXITCODE -ne 0) {
            throw "FreeFem++ failed with exit code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }

    $stepRows = @(Import-Csv -LiteralPath $stepData)
    if ($stepRows.Count -ne $Iterations) {
        throw "Expected $Iterations step rows; found $($stepRows.Count)."
    }
    $orderRows = @(Import-Csv -LiteralPath $orderData)
    if ($orderRows.Count -lt 4) {
        throw 'Shadow-order output is unexpectedly empty.'
    }

    Write-Output "Saved $($stepRows.Count) single-step rows: $stepData"
    Write-Output "Saved $($orderRows.Count) shadow rows: $orderData"
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
