$ErrorActionPreference = "Stop"

$source = $PSScriptRoot
$repo = Join-Path $source "PZT_droplets_GitHub"

if (-not (Test-Path (Join-Path $repo ".git"))) {
    throw "Cannot find the GitHub checkout: $repo"
}

$remote = git -C $repo remote get-url origin
if ($LASTEXITCODE -ne 0) { throw "Cannot read the Git remote." }
if ($remote -notmatch '^https://github\.com/alireza7660/PZT_droplets(?:\.git)?/?$') {
    throw "Unexpected repository: $remote"
}

robocopy "$source" "$repo" /E /XD PZT_droplets_GitHub .git .venv venv .pyenv __pycache__ .ipynb_checkpoints /XF *.h5 *.hdf5 *.hdf *.b16 *.b16dat .env .gitignore /MAX:94371840 /R:1 /W:1
if ($LASTEXITCODE -ge 8) { throw "File copying failed." }

git -C $repo add .
if ($LASTEXITCODE -ne 0) { throw "Git staging failed." }

$changes = git -C $repo status --porcelain
if ($LASTEXITCODE -ne 0) { throw "Cannot read Git status." }

if ($changes) {
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    git -C $repo commit -m "Auto sync $timestamp"
    if ($LASTEXITCODE -ne 0) { throw "Git commit failed." }
}

git -C $repo push origin main
if ($LASTEXITCODE -ne 0) { throw "Git push failed. See the message above." }

Write-Host "PZT files uploaded successfully." -ForegroundColor Green
