NAME                           = "MeshMap"
VERSION                        = 1
CONTRIBUTORS                   = ["ellis.langford.19@ucl.ac.uk"]
LAST_MOD_DATE                  = "30.10.2025"
VERBOSE                        = False

PARAMETERS = {
        "surface_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing surface .stl files"
        },
        "mesh_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing mesh .vtk files"
        },
        "global_info_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing global mesh info files"
        },
        "regional_info_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing regional mesh info files"
        },
        "surface_info_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing surface info files"
        },
        "dwi_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing diffusion files"
        },
        "cbf_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing cerebral blood flow files"
        },
        "props_fpath" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to an optional properties file containing"
                        " optional parameters."
        },
        "mesh_coarseness" : {
            "type"    : int,
            "default" : 20,
            "help"    : "Mesh coarseness from -50 (coarse) to +50 (fine)"
        },
        "zero_index_inputs" : {
            "type"    : bool,
            "default" : False,
            "help"    : "If True, input files are treated as 0-index, else 1-index. "
                        "Default is True."
        },
        "zero_index_outputs" : {
            "type"    : bool,
            "default" : False,
            "help"    : "If True, output files are saved as 0-index, else 1-index. "
                        "Default is True."
        },
        "incl_idx_col" : {
            "type"    : bool,
            "default" : True,
            "help"    : "If True, output files contain an index column. "
                        "Default is False."
        },
        "adjust_labels_dwi" : {
            "type"    : bool,
            "default" : True,
            "help"    : "If True, labels are updated based on FA. "
                        "Default is True."
        },
        "adjust_outer_labels" : {
            "type"    : bool,
            "default" : True,
            "help"    : "If True, labels on outer surface are adjusted "
                        "based on intersection with surface triangles. "
                        "Default is True."
        },
        "generate_cbf_map" : {
            "type"    : bool,
            "default" : True,
            "help"    : "If True, a scalar CBF map is produced. "
                        "Default is True."
        },
        "generate_fa_map" : {
            "type"    : bool,
            "default" : True,
            "help"    : "If True, a scalar FA map is produced. "
                        "Default is True."
        }
}