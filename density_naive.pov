#version 3.7;

global_settings { assumed_gamma 1.0 max_trace_level 12 }
background { color rgb <0.08, 0.08, 0.10> }

#declare DF3_PATH = "output/density_naive.df3"; // ou density_naive.df3

#declare XMAX = 10;
#declare YMAX = 10;
#declare ZMAX = 20;

// ---------------- CAMERA ----------------
camera {
  location <60, 18, -65>
  look_at  <0, 0, 0>
  angle 45
  right x*image_width/image_height
}


// ---------------- LIGHTS ----------------
// Une lumière blanche + deux jaunes
light_source { < 140,  120, -140> color rgb <1,1,1> }
light_source { <   0,   25,   10> color rgb <1.0,0.9,0.55> }
light_source { <   0,  -25,   10> color rgb <1.0,0.9,0.55> }

// ============================================================
// GRILLES PLUS LISIBLES (contraste + épaisseur + position)
// ============================================================
#declare PL = 24;         // demi-largeur plaques
#declare YPL = 12.5;      // plus proche du volume => plus visible
#declare THICK = 0.08;

#declare CELL  = 1.2;     // carreaux un peu plus serrés
#declare LINEW = 0.20;    // lignes plus épaisses => beaucoup plus visibles

#declare GridMask = function(x,z){
  max(
    select(LINEW - abs(sin(pi*x/CELL)), 0, 1),
    select(LINEW - abs(sin(pi*z/CELL)), 0, 1)
  )
};

#declare GridTexture =
texture{
  pigment{
    function{ GridMask(x,z) }
    color_map{
      [0.0 color rgb <0.92,0.92,0.92>] // carreaux clairs
      [1.0 color rgb <0.02,0.02,0.02>] // lignes quasi noires
    }
  }
  finish {
    ambient 0.35    // plus visible même si peu éclairé
    diffuse 0.65
    specular 0.25
    roughness 0.02
  }
}

// Plaque haut
box{
  <-PL,  YPL-THICK, -PL>, <PL, YPL+THICK, PL>
  texture { GridTexture }
  no_shadow
}

// Plaque bas
box{
  <-PL, -YPL-THICK, -PL>, <PL, -YPL+THICK, PL>
  texture { GridTexture }
  no_shadow
}

// ============================================================
// VOLUME DF3 (réglé pour laisser voir repères/grilles)
// ============================================================
box {
  <0,0,0>, <1,1,1>
  pigment { rgbt 1 }
  hollow

  interior {
    media {
      method 3
      intervals 1
      samples 70, 70

      // Un peu moins agressif que précédemment
      emission   0.85
      absorption 0.12
      scattering { 1, 0.65 }

      density {
        density_file df3 DF3_PATH interpolate 1

        // Transparence un peu augmentée pour laisser voir les axes
        color_map {
          [0.00 rgbt <0,0,0, 1.00>]
          [0.03 rgbt <0,0,0, 1.00>]

          [0.06 rgbt <0.05,0.25,1.00, 0.97>]
          [0.16 rgbt <0.10,0.45,1.00, 0.85>]

          [0.45 rgbt <0.65,0.85,1.00, 0.50>]
          [0.75 rgbt <1.00,1.00,1.00, 0.24>]
          [1.00 rgbt <1.00,1.00,1.00, 0.14>]
        }
      }
    }
  }

  translate <-0.5, -0.5, -0.5>
  scale <2*XMAX, 2*YMAX, 2*ZMAX>
  rotate <0, 12, 0>
  translate <4, 0, 0>
}

// ============================================================
// AXES TRÈS CLAIRS (épais + flèches + finish lumineux)
// ============================================================
#declare AX  = 18;
#declare RAD = 0.18;   // plus épais
#declare TIP = 0.55;   // taille flèche

#declare AxeFinish = finish { ambient 0.60 diffuse 0.40 specular 0.30 roughness 0.01 }

// X (orange)
union{
  cylinder { <-AX,0,0>, < AX,0,0>, RAD pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  cone     { < AX,0,0>, TIP, <AX+1.4,0,0>, 0 pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  cone     { <-AX,0,0>, TIP, <-AX-1.4,0,0>, 0 pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  no_shadow
}

// Y (orange)
union{
  cylinder { <0,-AX,0>, <0, AX,0>, RAD pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  cone     { <0, AX,0>, TIP, <0,AX+1.4,0>, 0 pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  cone     { <0,-AX,0>, TIP, <0,-AX-1.4,0>, 0 pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  no_shadow
}

// Z (orange)
union{
  cylinder { <0,0,-24>, <0,0,24>, RAD pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  cone     { <0,0, 24>, TIP, <0,0,25.8>, 0 pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  cone     { <0,0,-24>, TIP, <0,0,-25.8>, 0 pigment{ color rgb <1.0,0.70,0.10> } finish{ AxeFinish } }
  no_shadow
}