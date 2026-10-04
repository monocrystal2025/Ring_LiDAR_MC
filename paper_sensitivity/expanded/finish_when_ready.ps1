$ErrorActionPreference='Stop'
$expandedFolder='D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/expanded'
$deadline=(Get-Date).AddHours(8)
while ((Get-Date) -lt $deadline) {
    $campaignLog=Join-Path $expandedFolder 'campaign.log'
    if ((Test-Path -LiteralPath $campaignLog) -and (Select-String -LiteralPath $campaignLog -Pattern 'EXPANDED_CAMPAIGN_COMPLETE' -Quiet)) {
        $finishArguments='-wait -nodesktop -nosplash -logfile "D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/expanded/finish.log" -batch "addpath(''D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/expanded''); finish_expanded_robustness;"'
        $finishProcess=Start-Process -FilePath 'C:/Program Files/MATLAB/R2025a/bin/matlab.exe' -ArgumentList $finishArguments -WindowStyle Hidden -PassThru
        $finishProcess.Id | Set-Content -LiteralPath (Join-Path $expandedFolder 'finish_process_id.txt')
        $finishProcess.WaitForExit()
        $finishProcess.ExitCode | Set-Content -LiteralPath (Join-Path $expandedFolder 'finish_exit_code.txt')
        exit $finishProcess.ExitCode
    }
    Start-Sleep -Seconds 15
}
'No simulation restart; timed out awaiting campaign completion.' | Set-Content -LiteralPath (Join-Path $expandedFolder 'finish_timeout.txt')
exit 2
