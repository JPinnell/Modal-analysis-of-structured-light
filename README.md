# Modal-analysis-of-structured-light
MatLab codes to accompany Pinnell, J., Nape, I., Sephton, B., Cox, M., Rodriguez-Fajardo, V., and Forbes, A. "Modal analysis of structured light: a practical tutorial" (2020) 

The code herein allows the user to:
- Compute the detection cell mask required for optical mode projections using a camera (see FindDetCell.m).
- Generate digital holograms for encoding the basis functions on a phase-only device with the appropriate amplitude scaling (see GenBasisHolo.m)
- Simulate what the optical field should look like at the detection plane, for diagnostic purposes (see SimOpticalOverlap.m)
- Simulate various optical modal decompositions for different initial structured light beams and different bases (see SimModalDecomp.m)

Codes to generate Hermite-Gaussian (HG.m), Bessel-Gaussian (BG.m) and Laguerre-Gaussian (LG.m) structured light modes are also given for convenience and are required to run the simulation codes. A look up table for the inverse of the sinc function is also given (SincInv.mat) and is required for GenBasisHolo.m.

Please consider acknowledging our work if you have utilised the code herein.
