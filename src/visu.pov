#include "colors.inc"

camera {
  location <20, 20, -30>
  look_at <0, 0, 0>
}

light_source { <50, 50, -50> White }

background { Black }

interior {
  media {
    emission <1, 1, 1> * 2 // Glow intensity
    density {
      density_file df3 "density.df3"
      interpolate 1
      map_type 0
    }
  }
}

// Container box for the media
box {
  <0, 0, 0>, <1, 1, 1>
  pigment { rgbt 1 } // Transparent container
  hollow
  translate <-0.5, -0.5, -0.5> // Center it
  scale <20, 20, 40> // Scale to match physical dimensions (10+10, 10+10, 20+20)
}

