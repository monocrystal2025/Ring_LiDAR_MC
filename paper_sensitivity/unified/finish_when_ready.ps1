$ErrorActionPreference = 'Stop'
$unifiedRoot = 'D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/unified'
$simulationLog = Join-Path $unifiedRoot 'results/simulation.log'
$deadline = (Get-Date).AddHours(6)
while ((Get-Date) -lt $deadline) {
    $loggedComplete = (Test-Path -LiteralPath $simulationLog) -and
        (Select-String -LiteralPath $simulationLog -Pattern 'UNIFIED_SIMULATION_COMPLETE' -Quiet)
    $expectedWidths = @(5..100 | Where-Object { $_ % 5 -eq 0 }) + @(125,150)
    $allRawFilesReady = $true
    foreach ($caseNumber in 1..13) {
        foreach ($width in $expectedWidths) {
            $rawPath = Join-Path $unifiedRoot ('results/case_{0:D2}_width_{1:D3}.mat' -f $caseNumber,$width)
            if (-not (Test-Path -LiteralPath $rawPath)) { $allRawFilesReady = $false; break }
        }
        if (-not $allRawFilesReady) { break }
    }
    if ($loggedComplete -or $allRawFilesReady) {
        # Allow the last manifest write to finish. MATLAB then checks its
        # complete matrix; existence alone never passes numerical validation.
        Start-Sleep -Seconds 3
        $finishArguments = '-wait -nodesktop -nosplash -logfile "D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/unified/results/finish.log" -batch "addpath(''D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/unified''); finish_unified_robustness;"'
        $finishProcess = Start-Process -FilePath 'C:/Program Files/MATLAB/R2025a/bin/matlab.exe' -ArgumentList $finishArguments -WindowStyle Hidden -PassThru
        $finishProcess.Id | Set-Content -LiteralPath (Join-Path $unifiedRoot 'results/finish_process_id.txt')
        $finishProcess.WaitForExit()
        $finishProcess.ExitCode | Set-Content -LiteralPath (Join-Path $unifiedRoot 'results/finish_exit_code.txt')
        exit $finishProcess.ExitCode
    }
    Start-Sleep -Seconds 15
}
'Timed out waiting for simulation completion; no simulation was restarted.' | Set-Content -LiteralPath (Join-Path $unifiedRoot 'results/finish_wait_timeout.txt')
exit 2
