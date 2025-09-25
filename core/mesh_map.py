## Top-level wrapper
# Imports
import os
import sys
import shutil
import glob
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
        # Directories from input
        global_mesh_dir = self.get_parameter("global_mesh_dir")
        regional_mesh_dir = self.get_parameter("regional_mesh_dir")
        surface_dir = self.get_parameter("surface_dir")
        dwi_dir = self.get_parameter("dwi_dir")
        cbf_dir = self.get_parameter("cbf_dir")

        # Check directories
        for _dir in [global_mesh_dir, regional_mesh_dir,
                     surface_dir, dwi_dir, cbf_dir]:
            if _dir:
                if not os.path.isdir(_dir):
                    self.helpers.errors(f"Input directory does not exist {_dir}")
    
        # Helper for checking & copying files
        def copy_required(pattern, src_dir, dst_dir, file_name):
            matches = glob.glob(os.path.join(src_dir, pattern))
            if not matches:
                self.helpers.errors(f"Missing required file: {file_name} ({pattern}) in {src_dir}")
            if len(matches) > 1:
                self.helpers.errors(f"Multiple matches found for {file_name}: {matches}")
            os.makedirs(dst_dir, exist_ok=True)
            dst_path = os.path.join(dst_dir, os.path.basename(matches[0]))
            shutil.copy(matches[0], dst_path)
            return dst_path

        # Global mesh files
        global_dst = os.path.join(self.input_dir, "global_mesh")
        mesh_node_coords = copy_required("*coords*.txt", global_mesh_dir, global_dst, "Global node coordinates")
        mesh_tetra_indices = copy_required("*indices*.txt", global_mesh_dir, global_dst, "Global tetra indices")
        mesh_tetra_neighbours = copy_required("*neighbours*.txt", global_mesh_dir, global_dst, "Global tetra neighbours")
    
        # Regional mesh files
        regional_dst = os.path.join(self.input_dir, "regional_meshes")
        regions = ["cerebrum", "cerebellum", "cerebellumWM", "brainstem", "WM"]
    
        for region in regions:
            region_src = os.path.join(regional_mesh_dir, region)
            region_dst = os.path.join(regional_dst, region)
            copy_required(f"*{region}_L_node_coords*.txt", region_src, region_dst, f"{region} L node coords")
            copy_required(f"*{region}_L_tetra_indices*.txt", region_src, region_dst, f"{region} L tetra indices")
            copy_required(f"*{region}_R_node_coords*.txt", region_src, region_dst, f"{region} R node coords")
            copy_required(f"*{region}_R_tetra_indices*.txt", region_src, region_dst, f"{region} R tetra indices")
    
        # Surface files
        if surface_dir:
            surface_dst = os.path.join(self.input_dir, "surface_files")
            copy_required("*inner*coords*.txt", surface_dir, surface_dst, "Inner surface node coords")
            copy_required("*inner*indices*.txt", surface_dir, surface_dst, "Inner surface face indices")
            copy_required("*outer*coords*.txt", surface_dir, surface_dst, "Outer surface node coords")
            copy_required("*outer*indices*.txt", surface_dir, surface_dst, "Outer surface face indices")
    
        # DWI files
        if dwi_dir:
            dwi_dst = os.path.join(self.input_dir, "dwi_files")
            copy_required("*tensor*.nii*", dwi_dir, dwi_dst, "DTI tensor")
            copy_required("*L1*.nii*", dwi_dir, dwi_dst, "DTI L1")
            copy_required("*L2*.nii*", dwi_dir, dwi_dst, "DTI L2")
            copy_required("*L3*.nii*", dwi_dir, dwi_dst, "DTI L3")
            copy_required("*FA*.nii*", dwi_dir, dwi_dst, "DTI FA")
            copy_required("*MD*.nii*", dwi_dir, dwi_dst, "DTI MD")
    
        # CBF files
        if cbf_dir:
            cbf_dst = os.path.join(self.input_dir, "cbf_files")
            copy_required("*.nii*", cbf_dir, cbf_dst, "CBF map")
    
        # Set attributes
        self.global_mesh_node_coords = mesh_node_coords
        self.global_mesh_tetra_indices = mesh_tetra_indices
        self.global_mesh_tetra_neighbours = mesh_tetra_neighbours
        self.regional_mesh_dir = regional_dst
        self.surface_dir = surface_dst
        self.dwi_dir = dwi_dst
        self.cbf_dir = cbf_dst
    
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