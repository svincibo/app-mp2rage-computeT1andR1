[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.262-blue.svg)](https://doi.org/10.25663/brainlife.app.262)

# app-mp2rage-computeT1andR1
This app estimates T1map and R1map from MP2RAGE sequence. Before using this app to compute the T1 and R2 maps, you must first denoise the mp2rage uni image (see app-mp2rage-denoiseUNI).

Currently, only MP2RAGE data collected on Siemens MRI scanners is supported.

# Author
- Sophia Vinci-Booher (svincibo@iu.edu)

Code was adapted for brainlife.io; original version of this code can be found here: 
https://github.com/JosePMarques/MP2RAGE-related-scripts

# Project director
- Franco Pestilli (pestilli@utexas.edu)

# Funding
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-SMA-2004877](https://img.shields.io/badge/NSF_SMA-2004877-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=2004877)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://grantome.com/grant/NIH/R01-EB029272-01)
[![NIH-T32-HD007475](https://img.shields.io/badge/NIH_T32-HD007475-green.svg)](https://www.nichd.nih.gov/grants-contracts/training-careers/extramural/institutional)

# Citations.
We kindly ask that you cite the following articles when publishing papers and code using this code.

### brainlife.io platform Citations

Avesani, P., McPherson, B., Hayashi, S. et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. Sci Data 6, 69 (2019). https://doi.org/10.1038/s41597-019-0073-y

MIT Copyright (c) 2020 brainlife.io The University of Texas at Austin and Indiana University

### MP2RAGE Citations

Caan, M. W. A., Bazin, P.-L., Marques, J. P., de Hollander, G., Dumoulin, S. O., & van der Zwaag, W. (2019). MP2RAGEME: T1 , T2* , and QSM mapping in one sequence at 7 tesla. Human Brain Mapping, 40(6), 1786–1798.

Marques, J. P., & Gruetter, R. (2013). New Developments and Applications of the MP2RAGE Sequence-Focusing the Contrast and High Spatial Resolution R 1 Mapping. PloS One. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069294

Marques, J. P., & Gruetter, R. (2013). New developments and applications of the MP2RAGE sequence--focusing the contrast and high spatial resolution R1 mapping. PloS One, 8(7), e69294.

Marques, J. P., Kober, T., Krueger, G., van der Zwaag, W., Van de Moortele, P.-F., & Gruetter, R. (2010). MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage, 49(2), 1271–1281.

### Dependencies

This App only requires singularity to run (https://singularity.lbl.gov/install-request). If you don't have singularity, you will need to install following dependencies:

MATLAB: https://www.mathworks.com/products/matlab.html

