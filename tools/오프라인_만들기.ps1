# 수업용 오프라인 사본을 만든다. 서버도 인터넷도 없이 index.html 더블클릭으로 열린다.
#   powershell -ExecutionPolicy Bypass -File tools\오프라인_만들기.ps1
$ErrorActionPreference = 'Stop'
$root = Split-Path $PSScriptRoot -Parent
$name = "생명과학실험2_온라인실습"
$stage = Join-Path $env:TEMP $name
$zip = Join-Path (Split-Path $root -Parent) "${name}_$(Get-Date -Format yyyy-MM-dd).zip"

# 주차별 안내문을 guides.js 에 다시 내장한다 (file:// 은 fetch 가 막힌다)
python (Join-Path $PSScriptRoot 'build-guides.py')

Remove-Item -Recurse -Force $stage -ErrorAction SilentlyContinue
Copy-Item -Recurse (Join-Path $root 'docs') $stage
Copy-Item (Join-Path $PSScriptRoot '먼저_읽기.txt') $stage
Remove-Item $zip -ErrorAction SilentlyContinue
Compress-Archive -Path $stage -DestinationPath $zip -CompressionLevel Optimal
Remove-Item -Recurse -Force $stage
Write-Output ("{0}  {1:N0} KB" -f $zip, ((Get-Item $zip).Length / 1KB))
