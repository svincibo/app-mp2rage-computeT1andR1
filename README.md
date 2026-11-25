<h1 align="center">Compute T1 and R1 maps from MP2RAGE</h1>

<p align="center">
  Estimates quantitative T1 and R1 maps from a denoised Siemens MP2RAGE acquisition with a lookup-table method (compiled MATLAB, adapted from MP2RAGE-related-scripts).
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.514"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.514-blue.svg" alt="Run on brainlife.io"></a>
</p>

## At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Estimates quantitative T1 and R1 maps from a denoised Siemens MP2RAGE acquisition with a lookup-table method (compiled MATLAB, adapted from MP2RAGE-related-scripts).**

Reads the INV1, INV2 and UNI images with their JSON sidecars, builds the MP2RAGE signal lookup table from the sequence parameters and converts the UNI image into T1 and R1 (1/T1) maps; QA plots of the B1 sensitivity and value histograms are added to the product. Run app-mp2rage-denoiseUNI first. Only Siemens MP2RAGE data are supported; slicesperslab must be given when it is not in the image header.

| | |
|---|---|
| **Inputs** | `neuro/anat/mp2rage` · `neuro/mask` (optional) |
| **Outputs** | `neuro/anat/qmri` |
| **Container** | `brainlife/mcr:2020a` |
| **Resources** | 1 CPUs |
| **Run it** | [doi.org/10.25663/brainlife.app.514](https://doi.org/10.25663/brainlife.app.514) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `neuro/anat/mp2rage` | `neuro/anat/mp2rage` (denoised) | yes | Denoised MP2RAGE dataset: INV1 and INV2 magnitude images, the UNI image and their JSON sidecars (sequence parameters are read from them). |
| `neuro/mask` | `neuro/mask` | no | Optional brain mask applied to the maps. |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `output` | `neuro/anat/qmri` | T1map and R1 map |

<details>
<summary><b>Parameters</b> (1)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `slicesperslab` | number |  | Slices Per Slab should be available in the MRI sequence protocol. It is not always available in the image's header information. If you are sure that Slices Per Slab is available in the image's header information, you can leave this blank and we will get the value from there. Otherwise, please look up this value in the MRI sequence protocol and enter it here. |

</details>

<!-- brainlife-app:end -->

## Running the App

### On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.514](https://doi.org/10.25663/brainlife.app.514), choose `neuro/anat/mp2rage`, `neuro/mask` from your project, set the parameters and press **Execute**. The results are saved to your project as `neuro/anat/qmri`.

### Locally

You need [Singularity/Apptainer](https://apptainer.org/) (the App runs in `brainlife/mcr:2020a`), [`jq`](https://jqlang.github.io/jq/).

1. Clone this repository.
2. Create `config.json` with the paths to your inputs and the parameters:

   ```json
   {
       "mag_inv1": "/path/to/neuro/anat/mp2rage/mag_inv1",
       "mag_inv2": "/path/to/neuro/anat/mp2rage/mag_inv2",
       "unit1": "/path/to/neuro/anat/mp2rage/unit1",
       "json_inv1": "/path/to/neuro/anat/mp2rage/json_inv1",
       "json_inv2": "/path/to/neuro/anat/mp2rage/json_inv2",
       "json_unit1": "/path/to/neuro/anat/mp2rage/json_unit1",
       "mask": "/path/to/neuro/mask/mask",
       "slicesperslab": "<number>"
   }
   ```

3. Run `./main`.

## Troubleshooting

- **mp2rage_computeT1andR1.m didn't produce a T1.nii.gz file**: Check that the INV1, INV2 and UNI JSON sidecars contain the sequence parameters and set slicesperslab if it is missing.

## MP2RAGE Citations

Caan, M. W. A., Bazin, P.-L., Marques, J. P., de Hollander, G., Dumoulin, S. O., & van der Zwaag, W. (2019). MP2RAGEME: T1 , T2* , and QSM mapping in one sequence at 7 tesla. Human Brain Mapping, 40(6), 1786–1798.

Marques, J. P., & Gruetter, R. (2013). New Developments and Applications of the MP2RAGE Sequence-Focusing the Contrast and High Spatial Resolution R 1 Mapping. PloS One. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069294

Marques, J. P., & Gruetter, R. (2013). New developments and applications of the MP2RAGE sequence--focusing the contrast and high spatial resolution R1 mapping. PloS One, 8(7), e69294.

Marques, J. P., Kober, T., Krueger, G., van der Zwaag, W., Van de Moortele, P.-F., & Gruetter, R. (2010). MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage, 49(2), 1271–1281.

## Citing

If you use this App, please cite it, the software it builds on and brainlife.io:

- Compute T1 and R1 maps from MP2RAGE. brainlife.io App. [doi:10.25663/brainlife.app.514](https://doi.org/10.25663/brainlife.app.514)
- Marques JP, Gruetter R. New developments and applications of the MP2RAGE sequence - focusing the contrast and high spatial resolution R1 mapping. PLoS One 8(7):e69294 (2013). [doi:10.1371/journal.pone.0069294](https://doi.org/10.1371/journal.pone.0069294)
- Caan MWA, et al. MP2RAGEME: T1, T2*, and QSM mapping in one sequence at 7 tesla. Hum Brain Map 40(6):1786-1798 (2019).
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## Authors

- Sophia Vinci-Booher

Code was adapted for brainlife.io; original version of this code can be found here: 
https://github.com/JosePMarques/MP2RAGE-related-scripts

**Project director**

- Franco Pestilli

## Funding

[![Wellcome-226486/Z/22/Z](https://img.shields.io/badge/Wellcome-226486%2FZ%2F22%2FZ-purple.svg)](https://wellcome.org/)
[![NIH-NINDS-U24NS140384](https://img.shields.io/badge/NIH_NINDS-U24NS140384-green.svg)](https://reporter.nih.gov/project-details/11261037)
[![NIH-NINDS-UM1NS132207](https://img.shields.io/badge/NIH_NINDS-UM1NS132207-green.svg)](https://reporter.nih.gov/project-details/11411328)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://reporter.nih.gov/project-details/10262925)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://reporter.nih.gov/project-details/10428625)
[![NIH-NIMH-R01MH126699](https://img.shields.io/badge/NIH_NIMH-R01MH126699-green.svg)](https://reporter.nih.gov/project-details/10460628)
[![NIH-NIDCD-R01DC020088](https://img.shields.io/badge/NIH_NIDCD-R01DC020088-green.svg)](https://reporter.nih.gov/project-details/11140284)
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-2203524](https://img.shields.io/badge/NSF_BCS-2203524-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2203524)
[![NSF-IIS-1636893](https://img.shields.io/badge/NSF_IIS-1636893-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-OAC-2148729](https://img.shields.io/badge/NSF_OAC-2148729-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148729)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-IIS-2148700](https://img.shields.io/badge/NSF_IIS-2148700-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148700)
[![NSF-OAC-1916518](https://img.shields.io/badge/NSF_OAC-1916518-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1916518)

<sub>brainlife.io, which runs this App, is publicly funded by Wellcome, the U.S. National Institutes of Health and the U.S. National Science Foundation; acknowledging that support helps keep it running.</sub>
