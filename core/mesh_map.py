# Imports
import os
import sys
import shutil
import glob
import meshio
import numpy as np
from core import Core

# Add the utils directory to the Python path
utils_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.append(utils_path)

# Import custom modules
try:
    from base_cog import BaseCog
    from utils import Utils
    from helpers import Helpers, Path

except Exception as e:
    print(f"Failed to import toolkit modules from /app/utils - {e}")
    sys.exit(1)

# Create a class for processing
class MeshMap(BaseCog):
    def __init__(self, **kwargs):
        """
        Instantiate the Solver class.
        """
        super().__init__(**kwargs)
        
        # Instantiate custom modules
        self.utils = Utils()
        self.helpers = Helpers()
        self.path = Path()

        # Load parameters from CLI or properties file
        self.load_parameters()

    def copy_inputs(self):
        """
        Validate required input files and copy them to the working directory
        """        
        # Check directories
        if not self.get_parameter("mesh_dir"):
            if not self.get_parameter("global_info_dir") or not self.get_parameter("regional_info_dir"):
                self.helpers.errors(f"If no mesh .vtk files are provided, global and regional"
                                    f" mesh information files must be provided")
            if self.parameters["adjust_outer_labels"] and not self.get_parameter("surface_info_dir"):
                self.helpers.errors(f"Surface information files or meshes must be provided if adjustment"
                                    f" of outer tetrahedral labels ('adjust_outer_labels') is required")

        
        for dirname in ["surface_dir", "mesh_dir", 
                        "global_info_dir", 
                        "regional_info_dir", 
                        "surface_info_dir", 
                        "dwi_dir", "cbf_dir"]:
            _dir = self.get_parameter(dirname)
            if _dir:
                if not os.path.isdir(_dir):
                    self.helpers.errors(f"Input directory does not exist {_dir}")

        ### Surface .stl Files ###
        if self.get_parameter("surface_dir"):
            self.surface_dir = os.path.join(self.input_dir, "surface_files")
            shutil.copytree(self.get_parameter("surface_dir"), self.surface_dir)
        
        ### Mesh Files ###
        if self.get_parameter("mesh_dir"):
            self.mesh_dir = os.path.join(self.input_dir, "meshes")
            shutil.copytree(self.get_parameter("mesh_dir"), self.mesh_dir)

        ### Global Mesh Info Files ###
        global_info_dir = self.get_parameter("global_info_dir")
        if global_info_dir:
            global_dst = os.path.join(self.input_dir, "global_mesh_info")
            os.makedirs(global_dst, exist_ok=True)
            shutil.copy(glob.glob(os.path.join(global_info_dir, "*coords*.txt"))[0], os.path.join(global_dst, "global_node_coords.txt"))
            shutil.copy(glob.glob(os.path.join(global_info_dir, "*indices*.txt"))[0], os.path.join(global_dst, "global_tetra_indices.txt"))
            shutil.copy(glob.glob(os.path.join(global_info_dir, "*neighbours*.txt"))[0], os.path.join(global_dst, "global_tetra_neighbours.txt"))
        
        ### Regional Mesh Info Files ###
        regional_info_dir = self.get_parameter("regional_info_dir")
        if regional_info_dir:
            regions = ["cerebrum", "cerebellum", "cerebellumWM", "brainstem", "cerebrumWM"]
            for region in regions:
                region_src = os.path.join(regional_info_dir, region)
                region_dst = os.path.join(self.input_dir, "regional_mesh_info", region)
                os.makedirs(region_dst, exist_ok=True)
                shutil.copy(glob.glob(os.path.join(region_src, "*L*coords*.txt"))[0], os.path.join(region_dst, f"{region}_L_node_coords.txt"))
                shutil.copy(glob.glob(os.path.join(region_src, "*L*indices*.txt"))[0], os.path.join(region_dst, f"{region}_L_tetra_indices.txt"))
                shutil.copy(glob.glob(os.path.join(region_src, "*R*coords*.txt"))[0], os.path.join(region_dst, f"{region}_R_node_coords.txt"))
                shutil.copy(glob.glob(os.path.join(region_src, "*R*indices*.txt"))[0], os.path.join(region_dst, f"{region}_R_tetra_indices.txt"))
        
        ### Surface Info Files ###
        surface_info_dir = self.get_parameter("surface_info_dir")
        if surface_info_dir:
            surface_dst = os.path.join(self.input_dir, "surface_info")
            os.makedirs(surface_dst, exist_ok=True)
            shutil.copy(glob.glob(os.path.join(surface_info_dir, "*inner*coords*.txt"))[0], os.path.join(surface_dst, f"inner_surface_node_coords.txt"))
            shutil.copy(glob.glob(os.path.join(surface_info_dir, "*inner*indices*.txt"))[0], os.path.join(surface_dst, f"inner_surface_face_indices.txt"))
            shutil.copy(glob.glob(os.path.join(surface_info_dir, "*outer*coords*.txt"))[0], os.path.join(surface_dst, f"outer_surface_node_coords.txt"))
            shutil.copy(glob.glob(os.path.join(surface_info_dir, "*outer*indices*.txt"))[0], os.path.join(surface_dst, f"outer_surface_face_indices.txt"))

        ### Diffusion Weighted Imaging ###
        dwi_dir = self.get_parameter("dwi_dir")
        # If DWI options True
        if self.get_parameter("adjust_labels_dwi") or self.get_parameter("generate_fa_map"):
            # Check DWI directory and copy inputs
            if dwi_dir:
                dwi_dst = os.path.join(self.input_dir, "dwi_files")
                os.makedirs(dwi_dst, exist_ok=True)
                shutil.copy(glob.glob(os.path.join(dwi_dir, "*tensor*.nii*"))[0], os.path.join(dwi_dst, f"dwi_tensor.nii.gz"))
                shutil.copy(glob.glob(os.path.join(dwi_dir, "*L1*.nii.gz"))[0], os.path.join(dwi_dst, f"dwi_L1.nii.gz"))
                shutil.copy(glob.glob(os.path.join(dwi_dir, "*L2*.nii.gz"))[0], os.path.join(dwi_dst, f"dwi_L2.nii.gz"))
                shutil.copy(glob.glob(os.path.join(dwi_dir, "*L3*.nii.gz"))[0], os.path.join(dwi_dst, f"dwi_L3.nii.gz"))
                shutil.copy(glob.glob(os.path.join(dwi_dir, "*FA*.nii.gz"))[0], os.path.join(dwi_dst, f"dwi_FA.nii.gz"))
                shutil.copy(glob.glob(os.path.join(dwi_dir, "*MD*.nii.gz"))[0], os.path.join(dwi_dst, f"dwi_MD.nii.gz"))
                self.dwi_dir = dwi_dst
            else:
                self.helpers.errors(f"DWI directory must be provided if adjust_labels_dwi "
                                    f"or generate_fa_map options are set to True")
    
        ### Cerebral Blood Flow Imaging ###
        cbf_dir = self.get_parameter("cbf_dir")
        # If CBF options True
        if self.get_parameter("generate_cbf_map"):
            # Check CBF directory and copy inputs
            if cbf_dir:
                cbf_dst = os.path.join(self.input_dir, "cbf_files")
                os.makedirs(cbf_dst, exist_ok=True)
                shutil.copy(glob.glob(os.path.join(cbf_dir, "*.nii.gz"))[0], os.path.join(cbf_dst, f"cbf_map.nii.gz"))
                self.cbf_dir = cbf_dst
            else:
                self.helpers.errors(f"CBF directory must be provided if generate_cbf_map"
                                    f" option is set to True")            
        return

    def core(self):
        """
        Core processing
        """
        self.helpers.plugin_log(f"Starting {self.config['NAME']} at {self.helpers.now_time()}")

        # Tidy up log files from previous runs
        self.helpers.tidy_up_logs()

        # Directories
        self.input_dir   = os.path.join(self.base_dir, "inputs")
        self.interim_dir = os.path.join(self.base_dir, "interim_outputs")
        self.log_dir     = os.path.join(self.base_dir, "logs")
        self.output_dir  = os.path.join(self.base_dir, "outputs")

        for _dir in [self.input_dir, self.interim_dir, 
                     self.log_dir, self.output_dir]:
            shutil.rmtree(_dir, ignore_errors=True)
            os.makedirs(_dir, exist_ok=True)

        # Record parameters
        self.helpers.log_options(self.parameters)

        # Copy inputs
        self.helpers.plugin_log("Copying inputs")
        self.copy_inputs()

        # Set environment
        self.plugin_env = os.environ.copy()

        # Core processing
        proc = Core(self)
        proc.run()

        # Complete
        self.helpers.plugin_log(f"{self.config['NAME']} completed at {self.helpers.now_time()}")
        self.helpers.log_success()

# Set off pipeline object
if __name__ == "__main__":
    # Init
    plugin = MeshMap()
    # Execute
    plugin.core()