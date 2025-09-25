NAME                           = "MeshMap"
VERSION                        = 1
CONTRIBUTORS                   = ["ellis.langford.19@ucl.ac.uk"]
LAST_MOD_DATE                  = "25.09.2025"
VERBOSE                        = False

PARAMETERS = {
        "global_mesh_dir" : {
            "type"    : str,
            "help"    : "Path to directory containing global mesh files"
        },
        "regional_mesh_dir" : {
            "type"    : str,
            "help"    : "Path to directory containing regional mesh files"
        },
        "surface_dir" : {
            "type"    : str,
            "default" : "",
            "help"    : "Path to directory containing surface files"
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
        "write_mesh_info" : {
            "type"    : bool,
            "default" : True,
            "help"    : "If True, a txt file containing mesh information is produced. "
                        "Default is True."
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