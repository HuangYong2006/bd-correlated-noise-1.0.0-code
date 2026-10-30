# Package tracked repository files; local working material is never included.
[CmdletBinding()]
param(
    [string] $Version = '1.0.0',
    [string] $OutDir = (Join-Path $PSScriptRoot 'dist'),
    [switch] $Full,
    [switch] $Force
)
$ErrorActionPreference = 'Stop'
$Root = $PSScriptRoot
$Kind = if ($Full) { 'fulldata' } else { 'code' }
$ZipPath = Join-Path $OutDir "bd-correlated-noise-$Version-$Kind.zip"
if ((Test-Path -LiteralPath $ZipPath) -and -not $Force) {
    throw "Archive exists: $ZipPath (use -Force)"
}
$paths = @(git -C $Root -c core.quotepath=false ls-files)
if ($LASTEXITCODE -ne 0) { throw 'Cannot read Git file inventory' }
$paths = @($paths | Where-Object {
    (Test-Path -LiteralPath (Join-Path $Root $_) -PathType Leaf) -and
    ($Full -or -not $_.StartsWith('data/raw/'))
})
if (-not $paths -or 'code/bd_core.py' -notin $paths) {
    throw 'Missing repository sources; refusing to create archive'
}
New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
Add-Type -AssemblyName System.IO.Compression
$stream = [IO.File]::Open($ZipPath, [IO.FileMode]::Create)
$archive = [IO.Compression.ZipArchive]::new($stream, [IO.Compression.ZipArchiveMode]::Create, $false)
try {
    foreach ($relative in $paths) {
        $entry = $archive.CreateEntry($relative.Replace('\', '/'), [IO.Compression.CompressionLevel]::Optimal)
        $outputStream = $entry.Open()
        $inputStream = [IO.File]::OpenRead((Join-Path $Root $relative))
        try { $inputStream.CopyTo($outputStream) }
        finally { $inputStream.Dispose(); $outputStream.Dispose() }
    }
} finally { $archive.Dispose(); $stream.Dispose() }
Write-Output "Wrote $($paths.Count) tracked files to $ZipPath"
