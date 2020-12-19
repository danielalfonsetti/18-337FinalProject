# 18337FinalProject

## Abstract

Global sensitivity analysis (GSA) seeks to understand how distributions of the input parameters relate to distributions of the output parameters, and is thus different from local sensitivity analysis which is concerned with measuring sensitivity of the output at specific points of the input. In this project we implement four global sensitivity methods in Julia that have not been implemented previously: fractional factorial sampling, RBD-FAST, EASI, and the delta moment-independent measure. Time and memory benchmarking are performed for each, and numerical experiments from the methods’ corresponding research papers are recreated to help assert correctness when applicable.

## References

[1] Saltelli, A. (2008). Global sensitivity analysis: The primer. Chichester: Wiley.

[2] Tarantola, S., D. Gatelli and T. Mara (2006). "Random Balance Designs for the Estimation of First Order Global Sensitivity Indices.” Reliability Engineering and System Safety, 91:6, 717-727.

[3] Plischke, E. (2010). "An effective algorithm for computing global sensitivity indices (EASI).” Reliability Engineering & System Safety, 95:4, 354-360.

[4] Plischke, E., E. Borgonovo, and C. L. Smith (2013). "Global sensitivity measures from given data." European Journal of Operational Research, 226:3, 536-550.

[5] Herman, J., W. Usher (2019). Sensitivity Analysis Library (SALib).
