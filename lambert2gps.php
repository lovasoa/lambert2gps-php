<?php
function lambert2gps($x, $y, $lambert) {
      $lamberts = array(
       "LambertI" => 0,
       "LambertII" => 1,
       "LambertIII" => 2,
       "LamberIV" => 3,
       "LambertIIExtend" => 4,
       "Lambert93" => 5
      );
      $index = $lamberts[$lambert];
      $ntabs =  array(0.7604059656, 0.7289686274, 0.6959127966, 0.6712679322, 0.7289686274, 0.7256077650);
      $ctabs =  array(11603796.98, 11745793.39, 11947992.52, 12136281.99, 11745793.39, 11754255.426);
      $Xstabs = array(600000.0, 600000.0, 600000.0, 234.358, 600000.0, 700000.0);
      $Ystabs = array(5657616.674, 6199695.768, 6791905.085, 7239161.542, 8199695.768, 12655612.050);

      $n  = $ntabs [$index];
      $c  = $ctabs [$index];            // En mètres
      $Xs = $Xstabs[$index];          // En mètres
      $Ys = $Ystabs[$index];          // En mètres
      $l0 = 0.0;                    //correspond à la longitude en radian de Paris (2°20'14.025" E) par rapport à Greenwich
      $e = 0.08248325676;           //e du NTF (on le change après pour passer en WGS)
      $eps = 0.00001;     // précision


      /***********************************************************
      *  coordonnées dans la projection de Lambert 2 à convertir *
      ************************************************************/
      $X = $x;
      $Y = $y;

      /*
       * Conversion Lambert 2 -> NTF géographique : ALG0004
       */
        $R = Sqrt((($X - $Xs) * ($X - $Xs)) + (($Y - $Ys) * ($Y - $Ys)));
        $g = Atan(($X - $Xs) / ($Ys - $Y));

        $l = $l0 + ($g / $n);
        $L = -(1 / $n) * Log(Abs($R / $c));


        $phi0 = 2 * Atan(Exp($L)) - (pi() / 2.0);
        $phiprec = $phi0;
        $phii = 2 * Atan((Pow(((1 + $e * Sin($phiprec)) / (1 - $e * Sin($phiprec))), $e / 2.0) * Exp($L))) - (pi() / 2.0);

        while (!(Abs($phii - $phiprec) < $eps)) {
                $phiprec = $phii;
                $phii = 2 * Atan((Pow(((1 + $e * Sin($phiprec)) / (1 - $e * Sin($phiprec))), $e / 2.0) * Exp($L))) - (pi() / 2.0);
        }

        $phi = $phii;

    /*
     * Conversion NTF géogra$phique -> NTF cartésien : ALG0009
     */
    $a = 6378249.2;
    $h = 100;         // En mètres

    $N = $a / (Pow((1 - ($e * $e) * (Sin($phi) * Sin($phi))), 0.5));
    $X_cart = ($N + $h) * Cos($phi) * Cos($l);
    $Y_cart = ($N + $h) * Cos($phi) * Sin($l);
    $Z_cart = (($N * (1 - ($e * $e))) + $h) * Sin($phi);

    /*
     * Conversion NTF cartésien -> WGS84 cartésien : ALG0013
     */

    // Il s'agit d'une simple translation
    $XWGS84 = $X_cart - 168;
    $YWGS84 = $Y_cart - 60;
    $ZWGS84 = $Z_cart + 320;


    /*
     * Conversion WGS84 cartésien -> WGS84 géogra$phique : ALG0012
     */
    
    $l840 = 0.04079234433;    // 0.04079234433 pour passer dans un référentiel par rapport au méridien
    // de Greenwich, sinon mettre 0
    
    $e = 0.08181919106;              // On change $e pour le mettre dans le système WGS84 au lieu de NTF
    $a = 6378137.0;

    $P = Sqrt(($XWGS84 * $XWGS84) + ($YWGS84 * $YWGS84));

    $l84 = $l840 + Atan($YWGS84 / $XWGS84);

    $phi840 = Atan($ZWGS84 / ($P * (1 - (($a * $e * $e))
                                / Sqrt(($XWGS84 * $XWGS84) + ($YWGS84 * $YWGS84) + ($ZWGS84 * $ZWGS84)))));

    $phi84prec = $phi840;

    $phi84i = Atan(($ZWGS84 / $P) / (1 - (($a * $e * $e * Cos($phi84prec))
            / ($P * Sqrt(1 - $e * $e * (Sin($phi84prec) * Sin($phi84prec)))))));

    while (!(Abs($phi84i - $phi84prec) < $eps))
    {
        $phi84prec = $phi84i;
        $phi84i = Atan(($ZWGS84 / $P) / (1 - (($a * $e * $e * Cos($phi84prec))
                / ($P * Sqrt(1 - (($e * $e) * (Sin($phi84prec) * Sin($phi84prec))))))));

    }

    $phi84 = $phi84i;

    return array($phi84 * 180.0 / pi(), $l84 * 180.0 / pi());
}

function lambertI ($x, $y) {return lambert2gps($x,$y,"LambertI");} 
function lambertII ($x, $y) {return lambert2gps($x,$y,"LambertII");} 
function lambertIII ($x, $y) {return lambert2gps($x,$y,"LambertIII");} 
function lamberIV ($x, $y) {return lambert2gps($x,$y,"LamberIV");} 
function lambertIIExtend ($x, $y) {return lambert2gps($x,$y,"LambertIIExtend");} 
function lambert93 ($x, $y) {return lambert2gps($x,$y,"Lambert93");} 
?>
