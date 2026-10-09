# Pack this repository for upload to a Git host and/or a data archive.
#
# Produces two archives in -OutDir:
#
#   bd-correlated-noise-<Version>-code.zip
#       Everything needed to run the code and reproduce the manuscript numbers:
#       sources, docs, figures, manuscript sources and the derived analysis
#       tables. Raw simulation archives (data/raw/**/*.npz) are EXCLUDED,
#       because they are hundreds of MB and fully regenerable.
#
#   bd-correlated-noise-<Version>-fulldata.zip   (only with -Full)
#       The same tree, but including every raw .npz cell. Use this for a Zenodo
#       deposit referenced by the manuscript's Data Availability statement.
#
# Usage:
#   .\pack_repo.ps1
#   .\pack_repo.ps1 -Full
#   .\pack_repo.ps1 -Full -Version 1.0.0 -OutDir D:\releases

[CmdletBinding()]
param(
    [string] $Version = '1.0.0',
    [string] $OutDir  = (Join-Path $PSScriptRoot 'dist'),
    [switch] $Full,
    [switch] $Force
)

$ErrorActionPreference = 'Stop'

$Root = $PSScriptRoot
$Name = 'bd-correlated-noise'
$Kind = if ($Full) { 'fulldata' } else { 'code' }
$ZipPath = Join-Path $OutDir "$Name-$Version-$Kind.zip"
$Stage = Join-Path ([System.IO.Path]::GetTempPath()) "$Name-$Version-$Kind-stage"

# ---------------------------------------------------------------- guard rails
if (Test-Path $ZipPath) {
    if (-not $Force) { throw "Refusing to overwrite existing archive: $ZipPath  (use -Force)" }
    Remove-Item $ZipPath -Force
}
if (Test-Path $Stage) { Remove-Item $Stage -Recurse -Force }

New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
New-Item -ItemType Directory -Path $Stage  -Force | Out-Null

Write-Host "Packing  : $Root"
Write-Host "Mode     : $Kind"
Write-Host "Archive  : $ZipPath"
Write-Host ''

# --------------------------------------------------------------- collect files
$files = Get-ChildItem -Path $Root -Recurse -File -Force | Where-Object {
    $rel = $_.FullName.Substring($Root.Length).TrimStart('\', '/')
    $rel = $rel -replace '\\', '/'

    # never ship version-control, caches, environments or local run artefacts
    if ($rel -like '.git/*')             { return $false }
    if ($rel -like '*/__pycache__/*')    { return $false }
    if ($rel -like '__pycache__/*')      { return $false }
    if ($rel -like '*.pyc' -or $rel -like '*.pyo') { return $false }
    if ($rel -like '.venv/*' -or $rel -like 'venv/*') { return $false }
    if ($rel -like '.dsh_tmp/*')         { return $false }
    if ($rel -like '.dsh_acl_recovery/*'){ return $false }
    if ($rel -like '.figpreview/*')      { return $false }
    if ($rel -like '.pytest_cache/*')    { return $false }
    if ($rel -like 'dist/*')             { return $false }
    if ($rel -like '.ruff_cache/*')      { return $false }
    if ($_.Name -like 'acl-report-*.jsonl'  ) { return $false }
    if ($_.Name -like 'acl-backup-*.json'   ) { return $false }
    if ($_.Name -like 'acl-backup-*.ps1'    ) { return $false }

    # in code-only mode, drop the regenerable raw simulation archives
    if (-not $Full) {
        if ($rel -like 'data/raw/*') { return $false }
    }
    return $true
}

if (-not $files) { throw "No files selected; refusing to build an empty archive." }

# --------------------------------------------------------------------- staging
$count = 0
$bytes = 0L
foreach ($f in $files) {
    $rel  = $f.FullName.Substring($Root.Length).TrimStart('\', '/')
    $dest = Join-Path $Stage $rel
    $dir  = Split-Path $dest -Parent
    if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }
    Copy-Item -LiteralPath $f.FullName -Destination $dest -Force
    $count++
    $bytes += $f.Length
}

Write-Host ("Staged   : {0} files, {1:N1} MB" -f $count, ($bytes / 1MB))

# ------------------------------------------------------------------- integrity
# A release that cannot import its own model is not a release.
$core = Join-Path $Stage 'code/bd_core.py'
if (-not (Test-Path $core)) { throw "code/bd_core.py missing from the staged tree." }
foreach ($must in @('README.md', 'LICENSE', 'CITATION.cff', 'requirements.txt')) {
    if (-not (Test-Path (Join-Path $Stage $must))) {
        Write-Warning "Expected release file not found in staging: $must"
    }
}

# ----------------------------------------------------------------------- write
# NOTE: Compress-Archive (and ZipFile::CreateFromDirectory) store paths with the
# platform separator, i.e. "docs\DATA.md". Windows tolerates that, but the ZIP
# specification requires "/" and Linux/macOS unzip, GitHub's web view and the CI
# runners do not: a backslash archive extracts as flat files with literal
# backslashes in their names. We therefore write the entries ourselves with
# forward slashes, which is correct on every platform.
Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem

function New-PortableZip {
    param(
        [Parameter(Mandatory)] [string]   $SourceDirectory,
        [Parameter(Mandatory)] [string]   $DestinationPath,
        [Parameter(Mandatory)] [string[]] $RelativePaths
    )

    if (Test-Path $DestinationPath) { Remove-Item $DestinationPath -Force }
    $fs  = [System.IO.File]::Open($DestinationPath, [System.IO.FileMode]::CreateNew)
    $zip = New-Object System.IO.Compression.ZipArchive(
               $fs, [System.IO.Compression.ZipArchiveMode]::Create, $false)
    try {
        foreach ($rel in $RelativePaths) {
            $src   = Join-Path $SourceDirectory $rel
            $entry = $rel -replace '\\', '/'
            $ze    = $zip.CreateEntry($entry, [System.IO.Compression.CompressionLevel]::Optimal)
            $out   = $ze.Open()
            try {
                $in = [System.IO.File]::OpenRead($src)
                try { $in.CopyTo($out) } finally { $in.Dispose() }
            } finally { $out.Dispose() }
        }
    } finally {
        $zip.Dispose()
        $fs.Dispose()
    }
}

$relativePaths = $files | ForEach-Object { $_.FullName.Substring($Root.Length).TrimStart('\', '/') }
New-PortableZip -SourceDirectory $Root -DestinationPath $ZipPath -RelativePaths $relativePaths

Remove-Item $Stage -Recurse -Force

$zip = Get-Item $ZipPath
Write-Host ''
Write-Host ("WROTE    : {0}" -f $zip.FullName)
Write-Host ("SIZE     : {0:N1} MB" -f ($zip.Length / 1MB))
if (-not $Full) {
    Write-Host ''
    Write-Host 'Note: raw data/raw/**.npz cells are excluded from this archive.'
    Write-Host 'Re-run with -Full to build the complete data archive for a Zenodo deposit.'
}
