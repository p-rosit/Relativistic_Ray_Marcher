# Relativistic Ray Marcher

This is a project I made at LTH as a final project in EDAN35 High-Performance Computer Graphics.

This project tries to be a physically accurate ray marcher which renders images in the presence of a black hole which bends light in the ambient space. The black hole is assumed to induce a metric which bends light according to the Schwazschild metric which is the metric induced by a non-rotating black hole.

The project has not been modified and the bad practices are have not been fixed. For example, to image the scene all relevant parameters can be found in params.c but to change any render an image or several images with new parameters one needs to recompile the file ray_tracer.c. The project has been compiled with gcc 9.4.0 with the following parameters

```
gcc -fopenmp -O3 ray_tracer.c -lm -lpng
```

Textures are loaded as ppm files but images are saved as png files. This only due to the fact that I did not know how to load png files when the texture loading was implemented.
