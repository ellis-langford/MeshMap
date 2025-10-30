# Imports
import os
import sys
import shutil
import glob
import nibabel as nib
import numpy as np
import meshio
from scipy.spatial import KDTree
from scipy.sparse import coo_matrix
from collections import defaultdict

class Core(object):
    """
    Technical Setup.
    """
    def __init__(self, plugin_obj):
        # Check all expected attributed are present
        to_inherit = ["config", "utils", "helpers", "path", 
                      "parameters", "base_dir", "input_dir", 
                      "interim_dir", "output_dir",
                      "log_dir", "plugin_env"]
        for attr in to_inherit:
            try:
                setattr(self, attr, getattr(plugin_obj, attr))
            except AttributeError as e:
                print(f"Attribute Error - {e}")
                sys.exit(1)
        
    def read_txt(self, path, dtype=float, index_file=False, remove_idx_col=True):
        """
        Read txt file and optionally convert to zero-index and remove index column.
        
        Parameters:
        ---
        path (str) : Path to .txt file to be read
        dtype (type) : Type to read in file contents as
        index_file (bool) : If True, file contents is treated as indices and converted to zero-index if required
        remove_idx_col (bool) : If True, index column (if present), is removed from the file

        Returns:
        ---
        arr (np.array) : file contents saved as an array object
        """
        # Load file
        arr = np.loadtxt(path, dtype=dtype)

        if arr.ndim == 1:
            arr = arr[:, None]
    
        # Remove index column
        if remove_idx_col:
            first_col = arr[:, 0].astype(int)
            if np.all(first_col == np.arange(len(first_col))) or np.all(first_col == np.arange(1, len(first_col)+1)):
                arr = arr[:, 1:]
    
        # Convert to zero-index
        if not self.parameters["zero_index_inputs"] and index_file:
            if np.issubdtype(arr.dtype, np.integer):
                arr = arr - 1
            else:
                int_cols = np.all(arr == arr.astype(int), axis=0)
                arr[:, int_cols] = arr[:, int_cols].astype(int) - 1
    
        return arr

    def save_txt(self, path, data, fmt="%d", index_file=False, add_index_col=True):
        """
        Save numpy array to text file with optional index adjustment and index column.
        
        Parameters:
        ---
        path (str): Output file path
        data (np.array): Data to save
        fmt (str): Format string for np.savetxt
        index_file (bool): If True, adjusts for 1-indexing if required
        add_index_col (bool): If True, adds an index column to the output
        """
        # Ensure data is np array
        arr = np.asarray(data)
    
        if arr.ndim == 1:
            arr = arr[:, None]
    
        # Convert to one-indexing if required
        if not self.parameters["zero_index_outputs"] and index_file:
            if np.issubdtype(arr.dtype, np.integer):
                arr = arr + 1
            else:
                int_cols = np.all(arr == arr.astype(int), axis=0)
                arr[:, int_cols] = arr[:, int_cols].astype(int) + 1
    
        # Add index column if required
        if self.parameters["incl_idx_col"] and add_index_col:
            start = 0 if self.parameters["zero_index_outputs"] else 1
            idxs = np.arange(start, len(arr) + start)
            arr = np.column_stack([idxs, arr])

        # Save
        np.savetxt(path, arr, fmt=fmt)

        return

    def extract_mesh_info(self, mesh_path, output_dir, region, cell_type="tetra"):
        """
        Extract information from mesh files to txt files for mesh mapping.
    
        Parameters:
        ---
        mesh_path (str) : Path to directory containing .vtk mesh file
        output_dir (str) : Path to directory to save output .txt files to
        region (str) : Region that mesh belongs to, required for filename saving
        cell_type (str): Define whether input is mesh (tetra) or surface (triangles)
        """
        # Replace any 'unsigned_char' dtypes with 'int'
        if cell_type == "tetra":
            # Load file
            with open(mesh_path, "r") as f: 
                content = f.read()
    
            # Fix dtypes
            content = content.replace("unsigned_char", "int")
    
            # Save as txt
            base, _ = os.path.splitext(mesh_path)
            info_path = base + "_clean.vtk"
            with open(info_path, "w") as f: 
                f.write(content)

            mesh_path = info_path
        
        # Load mesh
        mesh = meshio.read(mesh_path)
        
        # Create output directory
        if region not in ["global", "outer_surface", "inner_surface"]:
            output_dir = os.path.join(output_dir, region)
        os.makedirs(output_dir, exist_ok=True)
    
        # Node coordinates and connectivity
        node_coords = mesh.points
        tetra = None
        element_type = ["tetra", "tetra10"] if cell_type == "tetra" else ["triangle", "tri"]
        for cellblock in mesh.cells:
            if cellblock.type in element_type:
                tetra = cellblock.data
                break

        unique_points, inverse = np.unique(tetra.flatten(), return_inverse=True)
        tetra = inverse.reshape(tetra.shape)
        node_coords = node_coords[unique_points]

        # Node coordinates file
        np.savetxt(os.path.join(output_dir, f"{region}_node_coords.txt"), node_coords, fmt="%.6f")

        # Node indices file
        if cell_type == "tetra":
            np.savetxt(os.path.join(output_dir, f"{region}_tetra_indices.txt"), tetra, fmt="%d")
        else:
            np.savetxt(os.path.join(output_dir, f"{region}_face_indices.txt"), tetra, fmt="%d")
    
        # Node neighbours file
        if region == "global":
            # Compute tetrahedron faces
            faces = np.vstack([
                np.sort(tetra[:, [0, 1, 2]], axis=1),
                np.sort(tetra[:, [0, 1, 3]], axis=1),
                np.sort(tetra[:, [0, 2, 3]], axis=1),
                np.sort(tetra[:, [1, 2, 3]], axis=1),
            ])
            # Tetrahedron IDs for each face
            tetra_ids = np.repeat(np.arange(len(tetra)), 4)
            
            # Sort vertex indices in each face
            faces_sorted = np.sort(faces, axis=1)
            
            # Build structured array for hashing and sorting
            dtype = np.dtype([('n0', faces_sorted.dtype),
                              ('n1', faces_sorted.dtype),
                              ('n2', faces_sorted.dtype)])
            faces_view = np.empty(faces_sorted.shape[0], dtype=dtype)
            faces_view['n0'] = faces_sorted[:, 0]
            faces_view['n1'] = faces_sorted[:, 1]
            faces_view['n2'] = faces_sorted[:, 2]
            
            # Sort faces to group identical together
            order = np.argsort(faces_view, order=('n0', 'n1', 'n2'))
            faces_sorted = faces_sorted[order]
            tetra_ids_sorted = tetra_ids[order]
            
            # Identify shared faces
            dupe_mask = np.all(faces_sorted[1:] == faces_sorted[:-1], axis=1)
            shared_face_indices = np.nonzero(dupe_mask)[0]
            
            # Build tetrahedron neighbour pairs
            t1 = tetra_ids_sorted[shared_face_indices]
            t2 = tetra_ids_sorted[shared_face_indices + 1]
            
            # Build symmetric adjacency
            all_pairs = np.vstack([np.column_stack([t1, t2]),
                                   np.column_stack([t2, t1])])
            
            # Sort and deduplicate
            all_pairs = np.unique(all_pairs, axis=0)
            
            # Construct adjacency list
            neighbours = defaultdict(list)
            for a, b in all_pairs:
                neighbours[a].append(b)
            
            # Convert to array form, padding with -1 (for boundary faces)
            max_neigh = max(len(v) for v in neighbours.values())
            tetra_neighbours = -np.ones((len(tetra), max_neigh), dtype=int)
            for k, v in neighbours.items():
                tetra_neighbours[k, :len(v)] = v
            
            # Save tetrahedral neighbours
            np.savetxt(os.path.join(output_dir, f"{region}_tetra_neighbours.txt"), tetra_neighbours, fmt="%d")

        return

    def load_directories(self):
        """
        Load required directories for processing
        """
        self.surface_dir = os.path.join(self.input_dir, "surface_files") if self.parameters["surface_dir"] else None
        self.mesh_dir = os.path.join(self.input_dir, "meshes") if self.parameters["mesh_dir"] else None
        self.global_info_dir = os.path.join(self.input_dir, "global_mesh_info") if self.parameters["global_info_dir"] else None
        self.regional_info_dir = os.path.join(self.input_dir, "regional_mesh_info") if self.parameters["regional_info_dir"] else None
        self.surface_info_dir = os.path.join(self.input_dir, "surface_info") if self.parameters["surface_info_dir"] else None
        self.dwi_dir = os.path.join(self.input_dir, "dwi_files") if self.parameters["dwi_dir"] else None
        self.cbf_dir = os.path.join(self.input_dir, "cbf_files") if self.parameters["cbf_dir"] else None

    def generate_mesh(self, region, mesh_coarseness=20):
        """
        THIS FUNCTION IS IN DEVELOPMENT
        DO NOT USE
        Generate a mesh .vtk from a surface .stl file
    
        Parameters:
        ---
        region (str) : Name of region to create mesh file for
        mesh_coarseness (int): Coarseness of mesh from -50 (coarse) to +50 (fine)
        """
        # Set up references
        # doc = sw.App.GetDocument()
        # model = doc.GetActiveModel()
        
        # # Import STL file
        # doc.ImportSurfaceFromStlFile(os.path.join(self.surface_files, f"{region}.stl", False, 1.0, False)
    
        # # Fix surface
        # doc.GetActiveSurface().Fix(SurfaceFixingControlParameters(9.9999999999999995e-07, 1.0000000000000001e-09))
    
        # # FE model creation
        # doc.CreateFeModel(f"{region} Model")
        
        # # Add surface to model
        # doc.EnableObjectsMode()
        # doc.GetModelByName(f"{region} Model").AddSurface(doc.GetSurfaceByName(region))
    
        # # Set model parameters
        # doc.EnableModelsMode()
        # model.SetExportType(Model.VtkVolume) # Set as VTK export
        # model.SetCompoundCoarsenessOnPart(model.GetPartByName(region), mesh_coarseness)
    
        # # Mesh generation
        # doc.GenerateMesh()
    
        # # Export mesh
        # doc.ExportVtkVolume(os.path.join(self.mesh_dir, f"{region}.vtk"), False)
    
        # # Disable region surface file
        # doc.GetSurfaceByName(region).SetVisible(False)
        # doc.EnableObjectsMode()

    def prepare_mesh_info_inputs(self):
        """
        Prepares input files by processing either surface files, mesh files or txt files
        """
        # Extract global mesh information
        global_mesh = glob.glob(os.path.join(self.mesh_dir, "global", "*global*.vtk"))[0]
        global_dst = os.path.join(self.input_dir, "global_mesh_info")
        self.extract_mesh_info(global_mesh, global_dst, "global")
        self.global_info_dir = os.path.join(self.input_dir, "global_mesh_info")

        # Extract regional mesh information
        regions = ["cerebrum_L", "cerebrum_R", "cerebrumWM_L", "cerebrumWM_R",
                   "cerebellum_L", "cerebellum_R", "cerebellumWM_L", "cerebellumWM_R",
                   "brainstem_L", "brainstem_R"]
        for region in regions:
            region_mesh = glob.glob(os.path.join(self.mesh_dir, region, f"*{region}*.vtk"))[0]
            regional_dst = os.path.join(self.input_dir, "regional_mesh_info")
            self.extract_mesh_info(region_mesh, regional_dst, region)
        self.regional_info_dir = os.path.join(self.input_dir, "regional_mesh_info")

        # Extract surface stl information
        if self.parameters["surface_dir"] and not self.parameters["surface_info_dir"]:
            surface_dst = os.path.join(self.input_dir, "surface_info")
            outer_surface = glob.glob(os.path.join(self.surface_dir, "*wholebrain*.stl"))[0]
            inner_surface = glob.glob(os.path.join(self.surface_dir, "*ventricles*.stl"))[0]
            self.extract_mesh_info(outer_surface, surface_dst, "outer_surface", cell_type="triangle")
            self.extract_mesh_info(inner_surface, surface_dst, "inner_surface", cell_type="triangle")
            self.surface_info_dir = os.path.join(self.input_dir, "surface_info")

        return
        
    def load_global_files(self):
        """
        Load global mesh coordinate, tetrahedra index and neighbour files.
        """   
        self.global_mesh_node_coords = self.read_txt(os.path.join(self.global_info_dir, "global_node_coords.txt"), dtype=float) # Node coords
        self.global_mesh_tetra_indices = self.read_txt(os.path.join(self.global_info_dir, "global_tetra_indices.txt"), dtype=int, index_file=True) # Tetrahedra node indices
        self.global_mesh_tetra_neighbours = self.read_txt(os.path.join(self.global_info_dir, "global_tetra_neighbours.txt"), dtype=int, index_file=True) # Shared tetrahedra faces

        return

    def match_faces(self):
        """
        Map inner and outer surface face nodes to global mesh node indices.
        """
        # Load surface files
        self.out_surface_nodes = self.read_txt(os.path.join(self.surface_info_dir, "outer_surface_node_coords.txt"), dtype=float)
        self.out_surface_faces = self.read_txt(os.path.join(self.surface_info_dir, "outer_surface_face_indices.txt"), dtype=int, index_file=True)
        
        # Build KDTree for nearest neighbour matching
        tree = KDTree(self.global_mesh_node_coords)
        dists, matched = tree.query(self.out_surface_nodes)

        # Check node matching against tolerances
        max_tol = 1e-3
        if np.any(dists > max_tol):
            self.helpers.plugin_log(f"WARNING: Some outer surface nodes mapped > {max_tol}. Max dist = {dists.max():.4e}")

        # Remap face indices to global indices
        remapped = matched[self.out_surface_faces]

        # Save output file
        output_file = os.path.join(self.interim_dir, f"outer_surface_face_indices_mapped.txt")
        self.save_txt(output_file, remapped, index_file=True, add_index_col=False)

        # Check file produced
        if not os.path.isfile(output_file):
            self.helpers.errors(f"Mapping surface face file not produced - {output_file}")

        return

    def classify_tetrahedra(self):
        """
        Classify which region each tetrahedral element belongs
        to based on circumsphere and centres
        """
        label_arrays = []
        region_files = []
    
        # Define regions
        regions = {
            "cerebrum": [1, 2],
            "cerebrumWM": [3, 4],
            "brainstem": [5, 6],
            "cerebellum": [7, 8],
            "cerebellumWM": [9, 10],
        }

        # Load node tree
        tree_global = KDTree(self.global_mesh_node_coords)

        # Loop over each region
        for region in regions:
            # Load region-specific coords and indices
            lh_node_coords   = self.read_txt(os.path.join(self.regional_info_dir, f"{region}_L", f"{region}_L_node_coords.txt"), dtype=float)
            lh_tetra_indices = self.read_txt(os.path.join(self.regional_info_dir, f"{region}_L", f"{region}_L_tetra_indices.txt"), dtype=int, index_file=True)
            rh_node_coords   = self.read_txt(os.path.join(self.regional_info_dir, f"{region}_R", f"{region}_R_node_coords.txt"), dtype=float)
            rh_tetra_indices = self.read_txt(os.path.join(self.regional_info_dir, f"{region}_R", f"{region}_R_tetra_indices.txt"), dtype=int, index_file=True)

            l_dists, _ = tree_global.query(lh_node_coords, k=1)
            r_dists, _ = tree_global.query(rh_node_coords, k=1)
            tol_L, tol_R = float(np.percentile(l_dists, 99)), float(np.percentile(r_dists, 99))

            # Build KD-trees
            tree_L = KDTree(lh_node_coords)
            tree_R = KDTree(rh_node_coords)

            labels = np.zeros(len(self.global_mesh_node_coords), dtype=np.int8)

            # Process in chunks to manage memory
            chunk_size = 500_000
            for start in range(0, len(self.global_mesh_node_coords), chunk_size):
                end = min(start + chunk_size, len(self.global_mesh_node_coords))
                chunk = self.global_mesh_node_coords[start:end]
        
                # Query distances
                dist_L, _ = tree_L.query(chunk, k=1)
                dist_R, _ = tree_R.query(chunk, k=1)
        
                # Assign labels
                labels_chunk = np.zeros(len(chunk), dtype=np.int8)
                labels_chunk[dist_L < tol_L] = regions[region][0]
                labels_chunk[dist_R < tol_R] = regions[region][1]
        
                labels[start:end] = labels_chunk
        
            # Save output
            output_file = os.path.join(self.interim_dir, f"{region}_labels.txt")
            self.save_txt(output_file, labels, add_index_col=False)

            label_arrays.append(labels)
            region_files.append(output_file)

        # Combine by taking the maximum label at each node
        n_nodes = len(label_arrays[0])
        combined_labels = np.zeros(n_nodes, dtype=int)
        
        for region, labels in zip(regions.keys(), label_arrays):
            mask = labels != 0
            combined_labels[mask] = labels[mask]

        # Identify labeled and unlabeled nodes
        unlabeled_idx = np.where(combined_labels == 0)[0]
        labeled_idx = np.where(combined_labels != 0)[0]
    
        # Build KD-tree on labeled nodes
        tree = KDTree(self.global_mesh_node_coords[labeled_idx])
    
        # Find nearest labeled node for each unlabeled one
        _, nearest = tree.query(self.global_mesh_node_coords[unlabeled_idx], k=1)
    
        # Assign the same label
        combined_labels[unlabeled_idx] = combined_labels[labeled_idx[nearest]]

        # Save combined labels
        output_file = os.path.join(self.interim_dir, "regional_labels.txt")
        self.save_txt(output_file, combined_labels, add_index_col=False)
        self.labels_file = output_file

        log_file = os.path.join(self.log_dir, "labelled_node_counts.txt")
        all_labels = np.unique(combined_labels)
        unlabelled = np.sum(combined_labels == 0)
        with open(log_file, 'w') as f:
            f.write(f"Total global nodes: {n_nodes:,}\n")
            f.write(f"Unlabelled nodes:   {unlabelled:,} ({unlabelled / n_nodes:.2%})\n")
    
        for region, (left_idx, right_idx) in regions.items():
            n_left = np.sum(combined_labels == left_idx)
            n_right = np.sum(combined_labels == right_idx)
            total = n_left + n_right
            with open(log_file, 'a') as f:
                f.write(f"{region:20s} L={n_left:,}  R={n_right:,}  Total={total:,}\n")

        # Check output file produced
        if not os.path.isfile(output_file):
            self.helpers.errors(f"Tetrahedra region label file not produced - {output_file}")

        return

    def revise_labels_by_dwi(self):
        """
        Revise mesh labels based on DWI FA values.
        """
        # Load DWI images (tensor components, eigenvalues, FA, MD)
        nii_files = {
            "tensor": os.path.join(self.dwi_dir, "dwi_tensor.nii.gz"),
            "L1": os.path.join(self.dwi_dir, "dwi_L1.nii.gz"),
            "L2": os.path.join(self.dwi_dir, "dwi_L2.nii.gz"),
            "L3": os.path.join(self.dwi_dir, "dwi_L3.nii.gz"),
            "FA": os.path.join(self.dwi_dir, "dwi_FA.nii.gz"),
            "MD": os.path.join(self.dwi_dir, "dwi_MD.nii.gz"),
        }
        imgs = {k: nib.load(v).get_fdata().astype(np.float32) for k, v in nii_files.items()}
    
        # Compute inverse affine for voxel mapping (world -> voxel space)
        ref_img = nib.load(nii_files["FA"])
        inv_affine = np.linalg.inv(ref_img.affine)
    
        # Compute centroids of all global tetrahedra
        coords = self.global_mesh_node_coords[self.global_mesh_tetra_indices, :]
        centroids = coords.mean(axis=1)
    
        # Map centroids to voxel coordinates and clip to bounds
        homog_centroids = np.c_[centroids, np.ones(len(centroids))]
        voxel_coords = (inv_affine @ homog_centroids.T).T[:, :3]
        voxel_coords = np.round(voxel_coords).astype(int)
    
        dimX, dimY, dimZ = imgs["FA"].shape
        voxel_coords[:, 0] = np.clip(voxel_coords[:, 0], 0, dimX - 1)
        voxel_coords[:, 1] = np.clip(voxel_coords[:, 1], 0, dimY - 1)
        voxel_coords[:, 2] = np.clip(voxel_coords[:, 2], 0, dimZ - 1)
    
        # Initialise array to store image values per tetrahedron
        N = self.global_mesh_tetra_indices.shape[0]
        tetra_center = np.zeros((N, 11), dtype=np.float32)

        # Sample DWI tensor components into first six columns
        for j in range(6):
            tetra_center[:, j] = imgs["tensor"][voxel_coords[:,0], voxel_coords[:,1], voxel_coords[:,2], j]

        # Sample L1, L2, L3, FA and MD into remaining columns
        tetra_center[:, 6] = imgs["L1"][voxel_coords[:,0], voxel_coords[:,1], voxel_coords[:,2]]
        tetra_center[:, 7] = imgs["L2"][voxel_coords[:,0], voxel_coords[:,1], voxel_coords[:,2]]
        tetra_center[:, 8] = imgs["L3"][voxel_coords[:,0], voxel_coords[:,1], voxel_coords[:,2]]
        tetra_center[:, 9] = imgs["FA"][voxel_coords[:,0], voxel_coords[:,1], voxel_coords[:,2]]
        tetra_center[:, 10] = imgs["MD"][voxel_coords[:,0], voxel_coords[:,1], voxel_coords[:,2]]
    
        # Clean up negative values and invalid FA ranges
        for col in [0, 3, 5, 6, 7, 8, 10]:
            tetra_center[tetra_center[:, col] < 0, col] = 0
        tetra_center[tetra_center[:, 9] > 1, 9] = 0
    
        # Neighbor averaging to fill gaps (using sparse neighbour matrix)
        rows = self.global_mesh_tetra_neighbours[:, 0]
        cols = self.global_mesh_tetra_neighbours[:, 5]
        data = np.ones_like(rows, dtype=np.float32)
        W = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

        # Compute neighbour counts and mean values
        neighbour_counts = np.array(W.sum(axis=1)).flatten()
        neighbour_sums = W @ tetra_center
        neighbour_means = np.divide(
            neighbour_sums,
            neighbour_counts[:, None],
            out=np.zeros_like(neighbour_sums),
            where=neighbour_counts[:, None] > 0,
        )

        # Replace zero entries with neighbour means
        mask = (tetra_center == 0)
        tetra_center[mask] = neighbour_means[mask]
    
        # Replace any remaining zeros with adjusted global mean
        for col in range(tetra_center.shape[1]):
            zeros = np.where(tetra_center[:, col] == 0)[0]
            if zeros.size:
                mean_val = tetra_center[:, col].mean()
                tetra_center[zeros, col] = mean_val + (mean_val * len(zeros)) / (N - len(zeros))
    
        # Save output arrays to files (tensor and FA)
        self.tensor_file = os.path.join(self.output_dir, "dwi_tensor.txt")
        self.FA_file = os.path.join(self.output_dir, "dwi_FA.txt")
        self.save_txt(self.tensor_file, tetra_center, fmt="%d " + "%f " * 11)
        self.save_txt(self.FA_file, tetra_center[:, 9].reshape(-1,1), fmt="%d %f")
    
        # Load labels for adjustment
        c = self.read_txt(self.tensor_file, dtype=np.float64, remove_idx_col=False)
        d = self.read_txt(self.labels_file, dtype=int, remove_idx_col=False)
    
        # Extract grey matter left (GML) and grey matter right (GMR) FA values
        GML = c[d[:,1] == 1][:, [0, 10]]
        GMR = c[d[:,1] == 2][:, [0, 10]]
        aveGML = np.mean(GML[:,1], dtype=np.float64) * 2.1
        aveGMR = np.mean(GMR[:,1], dtype=np.float64) * 2.1
    
        # Update labels based on FA thresholds
        for i in range(d.shape[0]):
            if d[i,1] in [1,3] and c[i,10] < aveGML:
                d[i,1] = 1
            if d[i,1] in [2,4] and c[i,10] < aveGMR:
                d[i,1] = 2
    
        # Save adjusted labels
        output_file = os.path.join(self.interim_dir, "labels_fa_adjusted.txt")
        self.save_txt(output_file, d, add_index_col=False)
        self.labels_file = output_file
    
        # Check output files were produced
        for fpath in [self.tensor_file, self.FA_file, output_file]:
            if not os.path.isfile(fpath):
                self.helpers.errors(f"File not produced - {fpath}")
    
        return

    def revise_outer_tetra_labels(self):
        """
        Revise tetrahedral labels for tetrahedra that are part of the outer surface.
        """
        # Load outer surface face indices and current labels
        outer_faces = self.read_txt(os.path.join(self.interim_dir, "outer_surface_face_indices_mapped.txt"), dtype=int, index_file=True)
        labels = self.read_txt(self.labels_file, dtype=int, remove_idx_col=False)
    
        # Sort vertices in each tetrahedron for consistent comparison
        tetra_sorted = np.sort(self.global_mesh_tetra_indices, axis=1)
        outer_sorted = np.sort(outer_faces, axis=1)
    
        # Add tetrahedron index as extra column for tracking
        tetra_sorted = np.hstack([tetra_sorted, np.arange(tetra_sorted.shape[0])[:, None]])
    
        # Sort rows lexicographically for fast intersection
        tetra_sorted = tetra_sorted[np.lexsort(tetra_sorted[:, ::-1].T)]
        outer_sorted = outer_sorted[np.lexsort(outer_sorted[:, ::-1].T)]
    
        # Helper: find intersecting rows by vertex combinations
        def intersect_rows(a, b):
            a = np.ascontiguousarray(a)
            b = np.ascontiguousarray(b)
            dtype = {"names": [f"f{i}" for i in range(a.shape[1])], "formats": a.shape[1]*[a.dtype]}
            return np.intersect1d(a.view(dtype), b.view(dtype), return_indices=True)[1]
    
        # Check all combinations of vertices
        rows_to_update = np.concatenate([
            intersect_rows(tetra_sorted[:, 0:3], outer_sorted), # First 3 vertices
            intersect_rows(tetra_sorted[:, 1:4], outer_sorted), # Last 3 vertices
            intersect_rows(tetra_sorted[:, [0,1,3]], outer_sorted), # Vertices 0,1,3
            intersect_rows(tetra_sorted[:, [0,2,3]], outer_sorted) # Vertices 0,2,3
        ])
    
        # Map back to original tetrahedron indices that need label updates
        tetra_indices = tetra_sorted[rows_to_update, 4].astype(int)
    
        # Update labels for outer face tetrahedra
        # 3&4 -> 1&2: cerebrum WM to cerebrum GM
        # 9&10 -> 7&8: cerebellum WM to cerebellum GM
        for old_label, new_label in zip([3, 4, 9, 10], [1, 2, 7, 8]):
            mask = labels[tetra_indices, 1] == old_label
            labels[tetra_indices[mask], 1] = new_label
    
        # Save revised labels
        output_file = os.path.join(self.interim_dir, "labels_outer_revised.txt")
        self.save_txt(output_file, labels, add_index_col=False)
    
        # Update attribute
        self.labels_file = output_file
    
        # Check output file produced
        if not os.path.isfile(output_file):
            self.helpers.errors(f"Revised outer label file not produced - {output_file}")
    
        return

    def map_scalar_to_tetra(self, nii_fpath, scalar_type):
        """
        Map a scalar field (e.g., CBF, FA) from a NIfTI volume onto the tetrahedral mesh.

        Parameters:
        ---
        nii_fpath (str): Path to the scalar NIfTI volume
        scalar_type (str): Name of the scalar type for output file naming (e.g., "CBF" or "FA")
        """
        # Load scalar NIfTI
        nii = nib.load(nii_fpath)
        img = nii.get_fdata()
        dimX, dimY, dimZ = img.shape # Dimensions of volume
        affine = nii.affine # Affine transformation matrix
        inv_affine = np.linalg.inv(affine) # Inverse for world -> coordinate mapping
    
        # Compute average non-zero voxels in image
        mask_nonzero = img != 0
        average_val = img[mask_nonzero].mean()
    
        # Compute centroids of each tetrahedron in global mesh coordinates
        idx = self.global_mesh_tetra_indices
        coords = self.global_mesh_node_coords[idx, :3]
        centroids = coords.mean(axis=1)
    
        # Convert centroids from world coordinates to voxel coordinates
        homog_centroids = np.c_[centroids, np.ones(len(centroids))]
        voxel_coords = (inv_affine @ homog_centroids.T).T[:, :3]
        voxel_coords = np.round(voxel_coords).astype(int)
    
        # Clip voxel coordinates to image bounds
        voxel_coords[:, 0] = np.clip(voxel_coords[:, 0], 0, dimX - 1)
        voxel_coords[:, 1] = np.clip(voxel_coords[:, 1], 0, dimY - 1)
        voxel_coords[:, 2] = np.clip(voxel_coords[:, 2], 0, dimZ - 1)
    
        # Sample image values at tetrahedral centroids
        tetra_vals = img[voxel_coords[:, 0], voxel_coords[:, 1], voxel_coords[:, 2]]
    
        # Fill zeros using neighborhood average
        zero_indices = np.where(tetra_vals == 0)[0]
        for i in zero_indices:
            x, y, z = voxel_coords[i]
            xs = np.clip(np.arange(x-2, x+1), 0, dimX-1)
            ys = np.clip(np.arange(y-2, y+1), 0, dimY-1)
            zs = np.clip(np.arange(z-2, z+1), 0, dimZ-1)
            cube = img[np.ix_(xs, ys, zs)]
            cube_flat = cube.flatten()
            mask_nonzero_cube = cube_flat != 0
            n_nonzero = np.sum(mask_nonzero_cube)
            n_zero = cube_flat.size - n_nonzero

            # If all values are non-zero, use global average
            if n_nonzero == 0:
                tetra_vals[i] = average_val
            # Else, use weighted mean of non-zero and average
            else:
                tetra_vals[i] = (cube_flat[mask_nonzero_cube].mean() * n_nonzero + n_zero * average_val) / cube_flat.size
    
        #Save scalar values
        output_file = os.path.join(self.output_dir, f"{scalar_type}_map.txt")
        self.save_txt(output_file, tetra_vals, fmt="%d %f")
    
        # Check output file produced
        if not os.path.isfile(output_file):
            self.helpers.errors(f"{scalar_type} scalar map file not produced - {output_file}")
    
        return

    def run(self):
        """
        Begin core processing.
        """
        self.helpers.plugin_log(f"Starting core processing at {self.helpers.now_time()}")

        # Load directories and attributes
        self.helpers.plugin_log(f"Loading directories")
        self.load_directories()

        regions = ["global", "ventricles", "brainstem_L", "brainstem_R",
                   "cerebrum_L", "cerebrum_R", "cerebrumWM_L", "cerebrumWM_R", 
                   "cerebellum_L", "cerebellum_R", "cerebellumWM_L", "cerebellumWM_R"]
        
        # If mesh info files not provided - generate
        if not self.parameters["global_info_dir"] and not self.parameters["regional_info_dir"]:
            self.helpers.plugin_log(f"Preparing mesh information inputs")
            self.prepare_mesh_info_inputs()

        # Read global mesh files
        self.helpers.plugin_log(f"Load global files")
        self.load_global_files()

        # Classify tetrahedra into regions
        self.helpers.plugin_log(f"Classifying tetrahedra into regions")
        self.classify_tetrahedra()
        self.labels_file = os.path.join(self.interim_dir, "regional_labels.txt")

        # Adjust labels based on FA
        if self.parameters["adjust_labels_dwi"]:
            self.helpers.plugin_log(f"Revising labels according to DWI FA values")
            self.revise_labels_by_dwi()

        # Revise outer tetrahedra labels
        if self.parameters["adjust_outer_labels"]:
            self.helpers.plugin_log(f"Adjusting outer tetrahedra labels")
            self.revise_outer_tetra_labels()

        # Produce CBF scalar map
        if self.parameters["generate_cbf_map"]:
            self.helpers.plugin_log(f"Generating CBF scalar map")
            cbf_file = glob.glob(os.path.join(self.cbf_dir, "*cbf_map*.nii*"))[0]
            self.map_scalar_to_tetra(cbf_file, "CBF")

        # Produce FA scalar map
        if self.parameters["generate_fa_map"]:
            self.helpers.plugin_log(f"Generating FA scalar map")
            fa_file = glob.glob(os.path.join(self.dwi_dir, "*dwi_FA*.nii*"))[0]
            self.map_scalar_to_tetra(fa_file, "FA")

        # Copy final label file
        shutil.copy(self.labels_file, os.path.join(self.output_dir, "labels.txt"))

        # Log completion
        self.helpers.plugin_log(f"Core analysis completed at {self.helpers.now_time()}")

        return