[![DOI](https://zenodo.org/badge/581099831.svg)](https://zenodo.org/badge/latestdoi/581099831)

# RBCFlowInRigidTube.m: Red blood cell motion and deformation inside of a rigid tube with varying cross-section.

- This code was developed for part of [my dissertation](https://www.researchgate.net/publication/355033649_Simulations_of_Red_Blood_Cell_Flow_by_Boundary_Integral_Methods) to simulate red blood cell flow using boundary integral methods.
- This repository contains the code for the examples presented in Chapter 6 for simulating the red blood cell motion and deformation in a microcapillary with constriction.

# Numerical examples

---

## Purely elastic membrane model

---

### Straight microcapillary

| Membrane shape | In-plane isotropic membrane tension |
| :-: | :-: |
|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/ElasRBC_Short_Pr4_2_Time0_75s/MembraneShapeInVesselElasRBC_Short_Pr4_2_Time0_75s_xy.gif">|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/ElasRBC_Short_Pr4_2_Time0_75s/isotropicTensionInVesselElasRBC_Short_Pr4_2_Time0_75s_xy.gif">|

|The fluid flow field inside and outside of the red blood cell|
| :-: |
|<video src="https://user-images.githubusercontent.com/13091572/218377243-10ea4ad7-82cd-4203-876b-6302f3812e9a.mp4" >|

---

### Short constricted microcapillary

| Membrane shape | In-plane isotropic membrane tension |
| :-: | :-: |
|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/ElasRBC_RefCons_6mic_Pr8/MembraneShapeInVesselElasRBC_RefCons_6mic_Pr8_xy.gif">|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/ElasRBC_RefCons_6mic_Pr8/isotropicTensionInVesselElasRBC_RefCons_6mic_Pr8_xy.gif">|

|The fluid flow field inside and outside of the red blood cell|
| :-: |
|<video src="https://user-images.githubusercontent.com/13091572/218379075-14bd8298-8632-4c65-a492-965457df0a0f.mp4">|
    
---
    
### Long constricted microcapillary

| Membrane shape | In-plane isotropic membrane tension |
| :-: | :-: |
|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/ElasRBC_LongConVes_Pr8/MembraneShapeInVesselElasRBC_LongConVes_Pr8_xy.gif">|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/ElasRBC_LongConVes_Pr8/isotropicTensionInVesselElasRBC_LongConVes_Pr8_xy.gif">|

|The fluid flow field inside and outside of the red blood cell|
| :-: |
|<video src="https://user-images.githubusercontent.com/13091572/218379238-23292720-f9b3-4854-b0cd-a9d047539c39.mp4">|
    
---

## Viscoelastic membrane model

---

### Straight microcapillary

| Membrane shape | In-plane isotropic membrane tension |
| :-: | :-: |
|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/MemVisRBC_Short_muMem10_Pr4_2_Time0_75/MembraneShapeInVesselMemVisRBC_Short_muMem10_Pr4_2_Time0_75_xy.gif">|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/MemVisRBC_Short_muMem10_Pr4_2_Time0_75/isotropicTensionInVesselMemVisRBC_Short_muMem10_Pr4_2_Time0_75_xy.gif">|

|The fluid flow field inside and outside of the red blood cell|
| :-: |
|<video src="https://user-images.githubusercontent.com/13091572/218379555-eaed5d43-3067-4a0f-9d6d-cfbd884f0145.mp4" >|


---

### Short constricted microcapillary

| Membrane shape | In-plane isotropic membrane tension |
| :-: | :-: |
|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/MemVisRBC_RefCons_6mic_muMem_3_18_Pr8/MembraneShapeInVesselMemVisRBC_RefCons_6mic_muMem_3_18_Pr8_xy.gif">|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/MemVisRBC_RefCons_6mic_muMem_3_18_Pr8/isotropicTensionInVesselMemVisRBC_RefCons_6mic_muMem_3_18_Pr8_xy.gif">|

|The fluid flow field inside and outside of the red blood cell|
| :-: |
|<video src="https://user-images.githubusercontent.com/13091572/218379713-2c17156b-61bc-4317-a42c-09d556f66d14.mp4">|

---
    
### Long constricted microcapillary

| Membrane shape | In-plane isotropic membrane tension |
| :-: | :-: |
|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/MemVisRBC_LongConVes_muMem_3_18_Pr40/MembraneShapeInVesselMemVisRBC_LongConVes_muMem_3_18_Pr40_xy.gif">|<img src="https://github.com/aligurbu/RBCFlowInRigidTube.m/blob/main/Results/MemVisRBC_LongConVes_muMem_3_18_Pr40/isotropicTensionInVesselMemVisRBC_LongConVes_muMem_3_18_Pr40_xy.gif">|

|The fluid flow field inside and outside of the red blood cell|
| :-: |
|<video src="https://user-images.githubusercontent.com/13091572/218379792-a80371ab-4aa1-4120-ad56-ab478b6d1b41.mp4">|

## Citation

    @phdthesis{gurbuz2021Thesis,
    title={Simulations of Red Blood Cell Flow by Boundary Integral Methods},
    author={G\"urb\"uz, Ali},
    year={2021},
    school={State University of New York at Buffalo}
    }
