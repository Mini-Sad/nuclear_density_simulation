#version 3.7;

// ============================================================
//  DENSITÉ NUCLÉAIRE 3D - RENDU PROFESSIONNEL
//  Projet INPS-2 - Visualisation volumétrique POV-Ray
// ============================================================

global_settings { 
  assumed_gamma 1.0 
  max_trace_level 15
  adc_bailout 0.01
  ambient_light rgb <0.3, 0.3, 0.35>
}

// Fond dégradé professionnel (noir vers bleu nuit)
sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 color rgb <0.02, 0.02, 0.06>]
      [0.3 color rgb <0.05, 0.05, 0.12>]
      [1.0 color rgb <0.08, 0.10, 0.18>]
    }
    scale 2
    translate -1
  }
}

#declare DF3_PATH = "output/density_optimized.df3";

#declare XMAX = 10;
#declare YMAX = 10;
#declare ZMAX = 20;

// ============================================================
//  CAMÉRA - Vue optimisée pour visualisation scientifique
// ============================================================
camera {
  location <55, 25, -60>
  look_at  <0, 0, 0>
  angle 42
  right x*image_width/image_height
  up y
}

// ============================================================
//  ÉCLAIRAGE - Setup studio professionnel à 3 points
// ============================================================
// Lumière principale (key light) - blanche douce
light_source { 
  <100, 80, -100> 
  color rgb <1.0, 1.0, 1.0> * 0.9
  area_light <15, 0, 0>, <0, 0, 15>, 5, 5
  adaptive 1
  jitter
}

// Lumière de remplissage (fill light) - bleutée
light_source { 
  <-80, 40, -60> 
  color rgb <0.6, 0.7, 0.9> * 0.5
  shadowless
}

// Lumière arrière (back/rim light) - chaude
light_source { 
  <0, 60, 80> 
  color rgb <1.0, 0.9, 0.7> * 0.4
  shadowless
}

// Lumière d'accentuation pour le volume
light_source { 
  <0, -30, 0> 
  color rgb <0.8, 0.85, 1.0> * 0.3
  shadowless
}

// ============================================================
//  GRILLES DE RÉFÉRENCE - Style scientifique élégant
// ============================================================
#declare PL = 22;
#declare YPL = 13;
#declare THICK = 0.05;
#declare CELL = 2.0;
#declare LINEW = 0.12;

// Fonction de masque pour grille
#declare GridMask = function(x,z){
  max(
    select(LINEW - abs(sin(pi*x/CELL)), 0, 1),
    select(LINEW - abs(sin(pi*z/CELL)), 0, 1)
  )
}

// Texture grille élégante semi-transparente
#declare GridTexture = texture {
  pigment {
    function { GridMask(x,z) }
    color_map {
      [0.0 color rgbt <0.85, 0.88, 0.95, 0.70>]  // fond semi-transparent
      [1.0 color rgbt <0.30, 0.35, 0.50, 0.40>]  // lignes plus opaques
    }
  }
  finish {
    ambient 0.4
    diffuse 0.5
    specular 0.15
    roughness 0.05
    reflection 0.02
  }
}

// Grille supérieure
box {
  <-PL, YPL-THICK, -PL>, <PL, YPL+THICK, PL>
  texture { GridTexture }
  no_shadow
}

// Grille inférieure
box {
  <-PL, -YPL-THICK, -PL>, <PL, -YPL+THICK, PL>
  texture { GridTexture }
  no_shadow
}

// ============================================================
//  VOLUME DF3 - Densité nucléaire avec colormap scientifique
// ============================================================
box {
  <0, 0, 0>, <1, 1, 1>
  pigment { rgbt 1 }
  hollow

  interior {
    media {
      method 3
      intervals 2
      samples 80, 80

      emission 0.75
      absorption 0.08
      scattering { 1, 0.45 }

      density {
        density_file df3 DF3_PATH interpolate 1

        // Colormap style "viridis" scientifique
        color_map {
          // Transparent pour les faibles densités
          [0.00 rgbt <0.267, 0.004, 0.329, 1.00>]  // violet foncé
          [0.02 rgbt <0.267, 0.004, 0.329, 0.98>]
          
          // Transition vers bleu
          [0.10 rgbt <0.282, 0.140, 0.458, 0.90>]  // violet-bleu
          [0.20 rgbt <0.254, 0.265, 0.530, 0.80>]  // bleu profond
          
          // Zone médiane cyan-vert
          [0.35 rgbt <0.191, 0.407, 0.556, 0.65>]  // bleu-cyan
          [0.50 rgbt <0.127, 0.566, 0.550, 0.50>]  // cyan-teal
          
          // Transition vers jaune-vert
          [0.65 rgbt <0.267, 0.679, 0.463, 0.35>]  // vert
          [0.80 rgbt <0.556, 0.759, 0.220, 0.25>]  // jaune-vert
          
          // Haute densité - jaune brillant
          [0.92 rgbt <0.850, 0.850, 0.098, 0.15>]  // jaune
          [1.00 rgbt <0.993, 0.906, 0.144, 0.08>]  // jaune vif
        }
      }
    }
  }

  translate <-0.5, -0.5, -0.5>
  scale <2*XMAX, 2*YMAX, 2*ZMAX>
  rotate <0, 8, 0>
  translate <2, 0, 0>
}

// ============================================================
//  AXES 3D - Style professionnel avec labels
// ============================================================
#declare AX = 18;
#declare RAD = 0.12;
#declare TIP = 0.40;

#declare AxeFinish = finish { 
  ambient 0.5 
  diffuse 0.4 
  specular 0.4 
  roughness 0.02
  metallic 0.3
}

// Couleurs des axes (RGB standard)
#declare ColorX = rgb <0.95, 0.25, 0.20>;  // Rouge (X)
#declare ColorY = rgb <0.20, 0.75, 0.25>;  // Vert (Y)  
#declare ColorZ = rgb <0.20, 0.45, 0.95>;  // Bleu (Z)

// Axe X (rouge)
union {
  cylinder { <-AX, 0, 0>, <AX, 0, 0>, RAD }
  cone { <AX, 0, 0>, TIP*1.2, <AX+1.2, 0, 0>, 0 }
  cone { <-AX, 0, 0>, TIP*1.2, <-AX-1.2, 0, 0>, 0 }
  pigment { color ColorX }
  finish { AxeFinish }
  no_shadow
}

// Axe Y (vert) - radial
union {
  cylinder { <0, -AX, 0>, <0, AX, 0>, RAD }
  cone { <0, AX, 0>, TIP*1.2, <0, AX+1.2, 0>, 0 }
  cone { <0, -AX, 0>, TIP*1.2, <0, -AX-1.2, 0>, 0 }
  pigment { color ColorY }
  finish { AxeFinish }
  no_shadow
}

// Axe Z (bleu) - axial
union {
  cylinder { <0, 0, -24>, <0, 0, 24>, RAD }
  cone { <0, 0, 24>, TIP*1.2, <0, 0, 25.5>, 0 }
  cone { <0, 0, -24>, TIP*1.2, <0, 0, -25.5>, 0 }
  pigment { color ColorZ }
  finish { AxeFinish }
  no_shadow
}

// ============================================================
//  MARQUEURS D'ÉCHELLE sur les axes
// ============================================================
#declare TickRad = 0.08;
#declare TickLen = 0.6;

// Ticks sur axe X
#for (i, -15, 15, 5)
  #if (i != 0)
    cylinder { 
      <i, -TickLen, 0>, <i, TickLen, 0>, TickRad 
      pigment { color ColorX * 0.8 }
      finish { AxeFinish }
      no_shadow
    }
  #end
#end

// Ticks sur axe Y
#for (i, -10, 10, 5)
  #if (i != 0)
    cylinder { 
      <-TickLen, i, 0>, <TickLen, i, 0>, TickRad 
      pigment { color ColorY * 0.8 }
      finish { AxeFinish }
      no_shadow
    }
  #end
#end

// Ticks sur axe Z
#for (i, -20, 20, 5)
  #if (i != 0)
    cylinder { 
      <0, -TickLen, i>, <0, TickLen, i>, TickRad 
      pigment { color ColorZ * 0.8 }
      finish { AxeFinish }
      no_shadow
    }
  #end
#end

// ============================================================
//  ORIGINE - Sphère centrale
// ============================================================
sphere {
  <0, 0, 0>, 0.35
  pigment { color rgb <1, 1, 1> }
  finish { 
    ambient 0.6 
    diffuse 0.3 
    specular 0.5 
    roughness 0.01
  }
  no_shadow
}

// ============================================================
//  PLAN DE COUPE SUBTIL (optionnel - visualisation symétrie)
// ============================================================
/*
disc {
  <0, 0, 0>, <1, 0, 0>, 15
  pigment { color rgbt <0.5, 0.6, 0.8, 0.85> }
  finish { ambient 0.3 diffuse 0.4 }
  no_shadow
}
*/