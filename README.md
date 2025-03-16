## MTU parameter identification
This code allows muscle model personalization by optimizing Muscle Tendon Unit (MTU) parameters.
The exeprimental setup is as follows:
* Joint (here ankle) torque from iskinetic dynamometer
        Filtered with ...
* Muscle activation from ElectroMyoGraphy (EMG)
        Rectified, Filtered with ..., normalized on Maximal Voluntary Contraction in Isometric
* Pennation angle and muscle fiber length from echography
        Extracted using custom code available [here](https://github.com/laboratoireIRISSE/EchoImageProcessing)
        
## Global overview
**Main.m** : Main code to run to optimize the MTU parameters (ℓmt, νmt; Fom, ℓom, ℓst, φo) from DeGroote's model. 
**ImprovedModel/main.m** : Main code to run to optimize the MTU parameters (...) from the new model we introduce in this study. 

## Code specifications
- MATLAB version: 9.13.0 (R2022b), Natick, Massachusetts: The MathWorks Inc.; 2022.
- The following code cannot run without "Casadi".