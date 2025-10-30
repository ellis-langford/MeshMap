<div align="center">
  <img src="./assets/mesh_map_logo.png" width="700">
  <br><br>
  <p align="center"><strong>An pipeline for mapping tissue classes onto tetrahedral meshes</strong></p>
</div>

<div align="center" style="display: flex; justify-content: center; gap: 10px; flex-wrap: wrap; margin-top: 10px;">
  <a href="https://profiles.ucl.ac.uk/101480-ellis-langford"><img src="https://custom-icon-badges.demolab.com/badge/UCL Profile-purple?logo=ucl" alt="UCL Profile"></a>
  <a href="https://orcid.org/0009-0006-1269-2632"><img src="https://img.shields.io/badge/ORCiD-green?logo=orcid&logoColor=white" alt="ORCiD"></a>
  <a href="https://github.com/ellis-langford"><img src="https://img.shields.io/badge/GitHub-%23121011.svg?logo=github&logoColor=white" alt="GitHub"></a>
  <a href="https://uk.linkedin.com/in/ellis-langford-8333441ab"><img src="https://custom-icon-badges.demolab.com/badge/LinkedIn-0A66C2?logo=linkedin-white&logoColor=fff" alt="LinkedIn"></a>
</div>

## Introduction

MeshMap is an image processing pipeline which determines which tissue class each element of a tetrahedral mesh belongs to.


## Requirements

To successfully run the MeshMap pipeline, please ensure the following requirements are met:

**Ubuntu 22.04 + Docker 27.3.1 + Python 3.10**<br>
*(other versions may be compatible but have not been tested)*


## Installation & Quick Start

To install the necessary components for MeshMap, please follow the steps below:

- Either, pull the docker image from GitHub container registry:

  ```bash
  docker pull ghcr.io/ellis-langford/mesh_map:v1
  ```

- Or clone the code from the GitHub repo and build image yourself:
  
  ```bash
  git clone https://github.com/ellis-langford/MeshMap.git
  cd MeshMap
  docker build -t ghcr.io/ellis-langford/mesh_map:v1 .
  ```
  
- Lauch a docker container from the MeshMap docker image:
  
  ```bash
  docker run -it -v /path/to/data:/path/to/data ghcr.io/ellis-langford/mesh_map:v1 bash
  ```

- Edit the example properties file to suit your requirements
  
  ```bash
  nano example_properties_file.json
  ```

- Navigate to your chosen output directory:
  
  ```bash
  cd /outputdir
  ```

- Run the pipeline:
  
  ```bash
  python3.10 /app/core/mesh_map.py --mesh_dir /path/to/mesh/dir --props_fpath /path/to/properties/file

## Inputs
1. `Global Mesh Files` & `Region Mesh Files`
   > A directory containing both global and regional mesh files can be supplied using *--mesh_dir*<br>

2. `Global Mesh Information Files` & `Region Mesh Information Files`
   > Can be used in place of a mesh directory - specified with *--global_info_dir* and *--regional_info_dir*<br>
   > *global_node_coords.txt:* List of node coordinates (x, y, z) for all mesh nodes.<br>
   > *global_tetra_indices.txt:* Tetrahedral connectivity: which 4 node indices form each tetrahedron.<br>
   > *global_tetra_neighbours.txt:* Neighbour tetrahedra for each element (used for interpolation).<br>
   > *{region}_L_node_coords.txt:* Node coordinates belonging to left-side ROI mesh.<br>
   > *{region}_L_tetra_indices.txt:* Connectivity (tetra indices) for the left ROI.<br>
   > *{region}_R_node_coords.txt:* Node coordinates belonging to right-side ROI mesh.<br>
   > *{region}_R_tetra_indices.txt:* Connectivity (tetra indices) for the right ROI.<br>

3. `Surface Files` (optional)
   > A directory surface .stl files or information files can be supplied using *--surface_dir* or *--surface_info_dir*<br>
   > *{region}.stl:* Surface .stl files for each region.<br>
   > *outer_surface_node_coords.txt:* Node coordinates of the outer brain surface.<br>
   > *outer_surface_face_indices.txt:* Surface triangle connectivity for outer surface.<br>

4. `Diffusion Weighted Imaging (DWI) Files` (optional)
   > A directory containing diffusion weighted imaging files can be supplied using *--dwi_dir*<br>
   > *dwi_tensor.nii.gz:* Diffusion tensor volume.<br>
   > *dwi_L1.nii.gz:* Principal eigenvalue.<br>
   > *dwi_L2.nii.gz:* Second eigenvalue.<br>
   > *dwi_L3.nii.gz:* Third eigenvalue.<br>
   > *dwi_FA.nii.gz:* Fractional anisotropy.<br>
   > *dwi_MD.nii.gz:* Mean diffusivity.<br>

5. `Cerebral Blood Flow (CBF) Files` (optional)
   > A directory containing cerebral blood flow files can be supplied using *--cbf_dir*<br>
   > *cbf_map.nii.gz:* CBF (cerebral blood flow) scalar field.<br>

## Pipeline Steps

After running these commands, the MeshMap pipeline will be run according to the options in the properties file. The pipeline steps available include:

1. `Match Faces`: Map outer surface nodes and faces to global indices.
2. `Classify Tetrahedra`: Tetrahedral elements are assigned a regional label.
3. `Revise Labels by DWI`: If flag *--adjust_labels_dwi* True, the labels will be updated according to DTI FA values.
4. `Revise Outer Tetra Labels`: If flag *--adjust_outer_labels* True, labels for tetrahedra on the outer surface are fixed.
5. `CBF Mapping`: If flag *--generate_cbf_map* True, a scalar CBF map is produced.
6. `FA Mapping`: If flag *--generate_fa_map* True, a scalar FA map is produced.

## Other Parameters
   > *zero_index_inputs:* If True, input files are expected to be zero-index<br>
   > *zero_index_outputs:* If True, output files are saved as zero-index<br>
   > *incl_idx_col:* If True, output files are saved containing an index column<br>

## Data Preparation

No additional data preparation is required.


## Output Structure

The output directory structure is as follows:

```
Output directory
├── inputs
│   ├── global_mesh
│   ├── regional_meshes
│   ├── surface_files
│   ├── dwi_files
│   └── cbf_files
├── interim_outputs
│   ├── outer_surface_face_indices_mapped.txt
│   ├── {region}_labels.txt
│   ├── regional_labels.txt
│   ├── labels_fa_adjusted.txt
│   └── labels_outer_revised.txt
├── logs
├── outputs
│   ├── labels.txt
│   ├── dwi_tensor.txt
│   ├── dwi_FA.txt
│   ├── CBF_map.txt
│   └── FA_map.txt
├── results.txt
└── errors.txt
```
- `inputs:` contains a copy of the input images
- `interim_outputs`: contains copies of files with various stages of processing applied
- `logs:` contains a plugin log (log.txt) and a record of the inputs and parameters (options.txt)
- `outputs:` contains the final output files
- `results.txt:` only produced if the pipeline executes successfully
- `errors.txt:` only produced if the pipeline fails to execute successfully (contains error info)


## Citation

```
@ARTICLE{xxxxxxxx,
  author={Langford E},
  journal={}, 
  title={}, 
  year={},
  volume={},
  number={},
  pages={},
  doi={}}
```
