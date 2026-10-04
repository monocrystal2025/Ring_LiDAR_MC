Set-StrictMode -Version Latest

function Sigmoid([double]$x) {
    return 1.0 / (1.0 + [Math]::Exp(-$x))
}

function New-SignalSequence([string]$Kind, [int]$N) {
    $S = New-Object double[] $N
    for ($i = 1; $i -le $N; $i++) {
        switch ($Kind) {
            'point' {
                $S[$i - 1] = 1.25 * [Math]::Exp(-0.5 * [Math]::Pow(($i - 25) / 5.0, 2))
            }
            'line' {
                $left = Sigmoid (($i - 14) / 1.6)
                $right = Sigmoid (($i - 36) / 1.6)
                $S[$i - 1] = 1.08 * ($left - $right)
            }
            'ring' {
                $lobe1 = 1.05 * [Math]::Exp(-0.5 * [Math]::Pow(($i - 16) / 3.2, 2))
                $lobe2 = 1.08 * [Math]::Exp(-0.5 * [Math]::Pow(($i - 34) / 3.4, 2))
                $S[$i - 1] = $lobe1 + $lobe2
            }
            default {
                throw "Unknown kind: $Kind"
            }
        }
    }
    return ,$S
}

function Add-Line([System.Collections.Generic.List[string]]$Svg, [double]$X1, [double]$Y1, [double]$X2, [double]$Y2, [string]$Style) {
    $Svg.Add(("<line x1='{0:F2}' y1='{1:F2}' x2='{2:F2}' y2='{3:F2}' style='{4}' />" -f $X1, $Y1, $X2, $Y2, $Style))
}

function Add-Rect([System.Collections.Generic.List[string]]$Svg, [double]$X, [double]$Y, [double]$W, [double]$H, [string]$Style) {
    $Svg.Add(("<rect x='{0:F2}' y='{1:F2}' width='{2:F2}' height='{3:F2}' style='{4}' />" -f $X, $Y, $W, $H, $Style))
}

function Add-Text([System.Collections.Generic.List[string]]$Svg, [double]$X, [double]$Y, [string]$Text, [int]$Size, [string]$Style) {
    $Escaped = [System.Security.SecurityElement]::Escape($Text)
    $Svg.Add(("<text x='{0:F2}' y='{1:F2}' font-size='{2}' style='{3}'>{4}</text>" -f $X, $Y, $Size, $Style, $Escaped))
}

function Add-Axis([System.Collections.Generic.List[string]]$Svg, [double]$Left, [double]$Base, [double]$Right, [double]$Top) {
    Add-Line $Svg $Left $Base $Right $Base "stroke:#09275a;stroke-width:2"
    Add-Line $Svg $Left $Base $Left $Top "stroke:#09275a;stroke-width:2"
    $Svg.Add(("<path d='M {0:F2} {1:F2} l -9 -5 l 0 10 z' fill='#09275a' />" -f $Right, $Base))
    $Svg.Add(("<path d='M {0:F2} {1:F2} l -5 9 l 10 0 z' fill='#09275a' />" -f $Left, $Top))
}

function New-AlignedPulseWindowSnrSvg([string]$Kind, [string]$OutFile) {
    $N = 49
    $M = 11
    $Threshold = 2.0
    $S = New-SignalSequence $Kind $N
    $Noise = New-Object double[] $N
    for ($i = 0; $i -lt $N; $i++) {
        $Noise[$i] = 0.16 + 0.015 * [Math]::Sin(($i + 1) * 0.65)
    }

    $W = New-Object double[] $N
    $Num = New-Object double[] $N
    $Den = New-Object double[] $N
    for ($i = 0; $i -lt $N; $i++) {
        $total = $S[$i] + $Noise[$i]
        if ($total -gt 0) {
            $W[$i] = $S[$i] / $total
            $Num[$i] = $W[$i] * $S[$i]
            $Den[$i] = $W[$i] * $W[$i] * $total
        }
    }

    $Dw = New-Object double[] $N
    $DecisionK = $N
    for ($k = 0; $k -lt $N; $k++) {
        $start = [Math]::Max(0, $k - $M + 1)
        $sumNum = 0.0
        $sumDen = 0.0
        for ($j = $start; $j -le $k; $j++) {
            $sumNum += $Num[$j]
            $sumDen += $Den[$j]
        }
        if ($sumDen -gt 0) {
            $Dw[$k] = $sumNum / [Math]::Sqrt($sumDen)
        }
        if (($DecisionK -eq $N) -and ($Dw[$k] -ge $Threshold)) {
            $DecisionK = $k
        }
    }
    if ($DecisionK -eq $N) {
        $DecisionK = [Array]::IndexOf($Dw, ($Dw | Measure-Object -Maximum).Maximum)
    }
    $WinStart = [Math]::Max(0, $DecisionK - $M + 1)

    $Width = 1700.0
    $Height = 980.0
    $Left = 170.0
    $Right = 1585.0
    $PlotW = $Right - $Left
    $BarW = 12.0
    $Rows = @(75.0, 290.0, 505.0, 720.0)
    $PlotH = 145.0
    $Svg = [System.Collections.Generic.List[string]]::new()
    $Svg.Add("<svg xmlns='http://www.w3.org/2000/svg' width='$Width' height='$Height' viewBox='0 0 $Width $Height'>")
    $Svg.Add("<rect width='100%' height='100%' fill='white' />")
    $Svg.Add("<style>text{font-family:Arial,Helvetica,sans-serif;fill:#071f4f}.small{font-size:20px}.label{font-weight:600}.formula{font-family:Times New Roman,serif;font-style:italic}</style>")

    function XPos([int]$Idx) {
        return $Left + (($Idx) * $PlotW / ($N - 1))
    }

    for ($r = 0; $r -lt 4; $r++) {
        Add-Rect $Svg 20 ($Rows[$r] - 45) 1660 195 "fill:none;stroke:#09275a;stroke-width:2;rx:14;ry:14"
        Add-Axis $Svg $Left ($Rows[$r] + $PlotH) $Right ($Rows[$r] + 5)
    }

    for ($i = 0; $i -lt $N; $i += 4) {
        $x = XPos $i
        Add-Line $Svg $x ($Rows[0] + 5) $x ($Rows[3] + $PlotH) "stroke:#b7c3d8;stroke-width:0.8;stroke-dasharray:4 7;opacity:0.45"
    }
    $xDecision = XPos $DecisionK
    Add-Line $Svg $xDecision ($Rows[0] + 5) $xDecision ($Rows[3] + $PlotH) "stroke:#c40046;stroke-width:1.6;stroke-dasharray:8 8;opacity:0.7"

    $titles = @{
        point = "Point-beam pulse sequence"
        line = "Line-beam pulse sequence"
        ring = "Ring-beam pulse sequence"
    }
    Add-Text $Svg 45 42 $titles[$Kind] 26 "font-weight:700"
    Add-Text $Svg 1240 42 "Same pulse index in all rows" 22 "font-weight:600;fill:#c40046"

    $rowLabels = @("1  Raw pulses", "2  Matched weight", "3  Same window terms", "4  Window statistic")
    for ($r = 0; $r -lt 4; $r++) {
        Add-Text $Svg 45 ($Rows[$r] + 10) $rowLabels[$r] 24 "font-weight:700"
    }

    $winX = XPos $WinStart
    $winRightX = XPos $DecisionK
    $winW = $winRightX - $winX
    for ($r = 0; $r -lt 3; $r++) {
        Add-Rect $Svg ($winX - $BarW) ($Rows[$r] + 12) ($winW + 2 * $BarW) ($PlotH - 12) "fill:#bfefff;stroke:#00a4df;stroke-width:2;opacity:0.28"
    }
    Add-Line $Svg $winX ($Rows[2] + 2) $winRightX ($Rows[2] + 2) "stroke:#09275a;stroke-width:2"
    Add-Line $Svg $winX ($Rows[2] + 2) ($winX + 12) ($Rows[2] - 5) "stroke:#09275a;stroke-width:2"
    Add-Line $Svg $winX ($Rows[2] + 2) ($winX + 12) ($Rows[2] + 9) "stroke:#09275a;stroke-width:2"
    Add-Line $Svg $winRightX ($Rows[2] + 2) ($winRightX - 12) ($Rows[2] - 5) "stroke:#09275a;stroke-width:2"
    Add-Line $Svg $winRightX ($Rows[2] + 2) ($winRightX - 12) ($Rows[2] + 9) "stroke:#09275a;stroke-width:2"
    Add-Text $Svg (($winX + $winRightX) / 2 - 18) ($Rows[2] - 12) "M" 25 "font-family:Times New Roman,serif;font-style:italic"

    $maxSN = 0.0
    for ($i = 0; $i -lt $N; $i++) {
        $maxSN = [Math]::Max($maxSN, $S[$i] + $Noise[$i])
    }
    $scale1 = 105.0 / $maxSN
    $base = $Rows[0] + $PlotH
    for ($i = 0; $i -lt $N; $i++) {
        $x = (XPos $i) - $BarW / 2
        $nH = $Noise[$i] * $scale1
        $sH = $S[$i] * $scale1
        Add-Rect $Svg $x ($base - $nH) $BarW $nH "fill:#b9bec7;stroke:#4a4f57;stroke-width:1"
        Add-Rect $Svg $x ($base - $nH - $sH) $BarW $sH "fill:#c40046;stroke:#8b0030;stroke-width:1"
    }

    $base = $Rows[1] + $PlotH
    $scale2 = 115.0
    for ($i = 0; $i -lt $N; $i++) {
        $x = (XPos $i) - $BarW / 2
        $h = $W[$i] * $scale2
        Add-Rect $Svg $x ($base - $h) $BarW $h "fill:#1f82d1;stroke:#0b3f7c;stroke-width:1"
    }

    $maxTerm = ($Num | Measure-Object -Maximum).Maximum
    $scale3 = 115.0 / $maxTerm
    $base = $Rows[2] + $PlotH
    for ($i = 0; $i -lt $N; $i++) {
        $inWin = ($i -ge $WinStart) -and ($i -le $DecisionK)
        $opacity = if ($inWin) { "1.0" } else { "0.22" }
        $x0 = XPos $i
        $hNum = $Num[$i] * $scale3
        $hDen = $Den[$i] * $scale3
        Add-Rect $Svg ($x0 - 6.5) ($base - $hDen) 5.5 $hDen "fill:#6585a8;stroke:#38506d;stroke-width:0.8;opacity:$opacity"
        Add-Rect $Svg ($x0 + 1.0) ($base - $hNum) 5.5 $hNum "fill:#c40046;stroke:#8b0030;stroke-width:0.8;opacity:$opacity"
    }

    $maxD = [Math]::Max(($Dw | Measure-Object -Maximum).Maximum, $Threshold * 1.2)
    $scale4 = 120.0 / $maxD
    $base = $Rows[3] + $PlotH
    $thresholdY = $base - $Threshold * $scale4
    Add-Line $Svg $Left $thresholdY $Right $thresholdY "stroke:#d00045;stroke-width:2;stroke-dasharray:10 8"
    Add-Text $Svg ($Left + 15) ($thresholdY - 10) "D_th = 2" 22 "fill:#d00045;font-family:Times New Roman,serif;font-style:italic"
    $path = ""
    for ($i = 0; $i -lt $N; $i++) {
        $x = XPos $i
        $y = $base - $Dw[$i] * $scale4
        if ($i -eq 0) {
            $path = "M {0:F2} {1:F2}" -f $x, $y
        } else {
            $path += (" L {0:F2} {1:F2}" -f $x, $y)
        }
    }
    $Svg.Add("<path d='$path' fill='none' stroke='#c40046' stroke-width='4' />")
    for ($i = 0; $i -lt $N; $i += 2) {
        $x = XPos $i
        $y = $base - $Dw[$i] * $scale4
        $Svg.Add(("<circle cx='{0:F2}' cy='{1:F2}' r='4.5' fill='#c40046' />" -f $x, $y))
    }
    $dy = $base - $Dw[$DecisionK] * $scale4
    $Svg.Add(("<circle cx='{0:F2}' cy='{1:F2}' r='8' fill='white' stroke='#c40046' stroke-width='4' />" -f $xDecision, $dy))
    Add-Text $Svg ($xDecision + 12) ($dy - 12) "first crossing" 20 "fill:#c40046;font-weight:600"

    Add-Text $Svg 230 ($Rows[0] + 18) "S_i: raw signal photons; N_i: noise photons" 22 ""
    Add-Text $Svg 230 ($Rows[1] + 18) "w_i = S_i / (S_i + N_i)" 24 "font-family:Times New Roman,serif;font-style:italic"
    Add-Text $Svg 230 ($Rows[2] + 18) "within the same W_k:  numerator w_i S_i;  variance term w_i^2(S_i+N_i)" 22 "font-family:Times New Roman,serif;font-style:italic"
    Add-Text $Svg 230 ($Rows[3] + 18) "D_w(k) = sum(w_i S_i) / sqrt(sum(w_i^2(S_i+N_i)))" 24 "font-family:Times New Roman,serif;font-style:italic"

    Add-Text $Svg ($Right - 115) ($Rows[0] + $PlotH + 35) "pulse index i" 21 ""
    Add-Text $Svg ($Right - 115) ($Rows[1] + $PlotH + 35) "pulse index i" 21 ""
    Add-Text $Svg ($Right - 115) ($Rows[2] + $PlotH + 35) "pulse index i" 21 ""
    Add-Text $Svg ($Right - 130) ($Rows[3] + $PlotH + 35) "window index k" 21 ""

    Add-Rect $Svg 1370 78 18 18 "fill:#c40046;stroke:#8b0030;stroke-width:1"
    Add-Text $Svg 1398 94 "S_i" 21 "font-family:Times New Roman,serif;font-style:italic"
    Add-Rect $Svg 1448 78 18 18 "fill:#b9bec7;stroke:#4a4f57;stroke-width:1"
    Add-Text $Svg 1476 94 "N_i" 21 "font-family:Times New Roman,serif;font-style:italic"
    Add-Rect $Svg 1370 522 18 18 "fill:#c40046;stroke:#8b0030;stroke-width:1"
    Add-Text $Svg 1398 538 "w_i S_i" 21 "font-family:Times New Roman,serif;font-style:italic"
    Add-Rect $Svg 1478 522 18 18 "fill:#6585a8;stroke:#38506d;stroke-width:1"
    Add-Text $Svg 1506 538 "w_i^2(S_i+N_i)" 21 "font-family:Times New Roman,serif;font-style:italic"

    $Svg.Add("</svg>")
    Set-Content -LiteralPath $OutFile -Value $Svg -Encoding UTF8
}

$AssetDir = Split-Path -Parent $MyInvocation.MyCommand.Path
New-AlignedPulseWindowSnrSvg "point" (Join-Path $AssetDir "aligned_point_pulse_window_snr_decision.svg")
New-AlignedPulseWindowSnrSvg "line" (Join-Path $AssetDir "aligned_line_pulse_window_snr_decision.svg")
New-AlignedPulseWindowSnrSvg "ring" (Join-Path $AssetDir "aligned_ring_pulse_window_snr_decision.svg")
