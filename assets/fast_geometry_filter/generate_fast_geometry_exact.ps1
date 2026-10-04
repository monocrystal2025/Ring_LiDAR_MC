Add-Type -AssemblyName System.Drawing

$outDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$W = 1600
$H = 900

function New-Figure {
    $bmp = [System.Drawing.Bitmap]::new($W, $H)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
    $g.Clear([System.Drawing.Color]::White)
    return @{ Bitmap = $bmp; Graphics = $g }
}

function Save-Figure($fig, $name) {
    $path = Join-Path $outDir $name
    $fig.Bitmap.Save($path, [System.Drawing.Imaging.ImageFormat]::Png)
    $fig.Graphics.Dispose()
    $fig.Bitmap.Dispose()
}

function New-Pen($color, $width) {
    return [System.Drawing.Pen]::new($color, $width)
}

function New-Brush($color) {
    return [System.Drawing.SolidBrush]::new($color)
}

function Draw-DashedCircle($g, $cx, $cy, $r) {
    $pen = New-Pen ([System.Drawing.Color]::FromArgb(255, 0, 92, 255)) 4
    $pen.DashPattern = [float[]](14, 11)
    $g.DrawEllipse($pen, $cx - $r, $cy - $r, 2 * $r, 2 * $r)
    $pen.Dispose()
}

function Draw-DashedRectangle($g, $x, $y, $w, $h) {
    $pen = New-Pen ([System.Drawing.Color]::FromArgb(255, 0, 92, 255)) 4
    $pen.DashPattern = [float[]](14, 11)
    $g.DrawRectangle($pen, $x, $y, $w, $h)
    $pen.Dispose()
}

function Get-LineSamples {
    param(
        [double]$x0,
        [double]$y0,
        [double]$x1,
        [double]$y1,
        [int]$n
    )
    $pts = New-Object System.Collections.Generic.List[object]
    for ($i = 0; $i -lt $n; $i++) {
        $t = [double]$i / [double]($n - 1)
        $px = $x0 + (($x1 - $x0) * $t)
        $py = $y0 + (($y1 - $y0) * $t)
        $pts.Add([double[]]@($px, $py)) | Out-Null
    }
    return $pts.ToArray()
}

function Draw-Dot($g, $x, $y, $accepted) {
    if ($accepted) {
        $brush = New-Brush ([System.Drawing.Color]::FromArgb(255, 0, 145, 32))
    } else {
        $brush = New-Brush ([System.Drawing.Color]::FromArgb(210, 105, 105, 105))
    }
    $r = 8
    $g.FillEllipse($brush, $x - $r, $y - $r, 2 * $r, 2 * $r)
    $brush.Dispose()
}

function Draw-PointFigure {
    $fig = New-Figure
    $g = $fig.Graphics
    $cx = 800
    $cy = 450
    $pointRadius = 335

    for ($i = 72; $i -ge 1; $i--) {
        $drawRadius = $pointRadius * $i / 72
        $u = $i / 72
        $alpha = [int](5 + 34 * [Math]::Pow(1 - $u, 1.6))
        $brush = New-Brush ([System.Drawing.Color]::FromArgb($alpha, 0, 102, 255))
        $g.FillEllipse($brush, $cx - $drawRadius, $cy - $drawRadius, 2 * $drawRadius, 2 * $drawRadius)
        $brush.Dispose()
    }

    Draw-DashedCircle $g $cx $cy $pointRadius
    $pts = Get-LineSamples 110 690 1490 230 33
    foreach ($p in $pts) {
        $px = [double]$p[0]
        $py = [double]$p[1]
        $d = [Math]::Sqrt([Math]::Pow($px - $cx, 2) + [Math]::Pow($py - $cy, 2))
        Draw-Dot $g $px $py ($d -le ($pointRadius - 2))
    }
    Save-Figure $fig 'fast_geometry_point_pulse_samples_exact.png'
}

function Draw-LineFigure {
    $fig = New-Figure
    $g = $fig.Graphics
    $x = 155
    $y = 295
    $w = 1290
    $h = 310
    $yc = $y + $h / 2

    for ($dy = $h / 2; $dy -ge 1; $dy -= 3) {
        $u = [Math]::Abs($dy) / ($h / 2)
        $alpha = [int](6 + 42 * [Math]::Exp(-3.2 * $u * $u))
        $brush = New-Brush ([System.Drawing.Color]::FromArgb($alpha, 0, 102, 255))
        $g.FillRectangle($brush, $x, $yc - $dy, $w, 2 * $dy)
        $brush.Dispose()
    }

        $centerPen = New-Pen ([System.Drawing.Color]::FromArgb(115, 255, 255, 255)) 7
    $g.DrawLine($centerPen, $x, $yc, $x + $w, $yc)
    $centerPen.Dispose()

    Draw-DashedRectangle $g $x $y $w $h
    $pts = Get-LineSamples 100 720 1500 185 35
    foreach ($p in $pts) {
        $px = [double]$p[0]
        $py = [double]$p[1]
        $accepted = ($px -ge $x) -and ($px -le ($x + $w)) -and ($py -ge $y) -and ($py -le ($y + $h))
        Draw-Dot $g $px $py $accepted
    }
    Save-Figure $fig 'fast_geometry_line_pulse_samples_exact.png'
}

function Draw-RingFigure {
    $fig = New-Figure
    $g = $fig.Graphics
    $cx = 800
    $cy = 450
    $inner = 190
    $outer = 395
    $r0 = ($inner + $outer) / 2
    $sigma = 52

    for ($r = $inner; $r -le $outer; $r += 4) {
        $q = ($r - $r0) / $sigma
        $alpha = [int](9 + 68 * [Math]::Exp(-0.5 * $q * $q))
        $pen = New-Pen ([System.Drawing.Color]::FromArgb($alpha, 0, 102, 255)) 8
        $g.DrawEllipse($pen, $cx - $r, $cy - $r, 2 * $r, 2 * $r)
        $pen.Dispose()
    }

    Draw-DashedCircle $g $cx $cy $inner
    Draw-DashedCircle $g $cx $cy $outer
    $pts = Get-LineSamples 115 725 1485 145 33
    foreach ($p in $pts) {
        $px = [double]$p[0]
        $py = [double]$p[1]
        $d = [Math]::Sqrt([Math]::Pow($px - $cx, 2) + [Math]::Pow($py - $cy, 2))
        Draw-Dot $g $px $py (($d -ge ($inner + 2)) -and ($d -le ($outer - 2)))
    }
    Save-Figure $fig 'fast_geometry_ring_pulse_samples_exact.png'
}

Draw-PointFigure
Draw-LineFigure
Draw-RingFigure
