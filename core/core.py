# Imports
import os
import sys
import shutil
import glob
import nibabel as nib
import numpy as np
from scipy.spatial import KDTree
from scipy.sparse import coo_matrix

class Core(object):
    """
    Technical Setup.
    """
    def __init__(self, plugin_obj):
        # Check all expected attributed are present
        to_inherit = ["config", "utils", "helpers", "path", "parameters",
                      "base_dir", "input_dir", "interim_dir", "output_dir",
                      "global_mesh_node_coords", "global_mesh_tetra_indices", 
                      "global_mesh_tetra_neighbours", "regional_mesh_dir", 
                      "surface_dir", "dwi_dir", "cbf_dir", "log_dir", "plugin_env"]
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

    def tetrahedron_circumscribed_sphere(self, T):
        """
        Calculate circumscribed sphere of a tetrahedron using least squares.
        
        Parameters:
        ---
        T (np.array): Coordinates of tetrahedron nodes (shape: 4x3)
        
        Returns:
        ---
        C (np.array): Circumsphere centre (3D coordinates)
        r (float): Circumsphere radius
        """
        A = 2 * (T[:, 1:] - T[:, [0]])
        b = np.sum(T[:, 1:] ** 2 - T[:, [0]] ** 2, axis=0)
        C = np.linalg.lstsq(A.T, b, rcond=None)[0]
        return C, np.linalg.norm(C - T[:, 0])

    def compute_tetra_centres(self, nodes, cell_labels):
        """
        Compute circumsphere centres and radii for a set of tetrahedra.
        
        Parameters:
        ---
        nodes (np.array): Node coordinates
        cell_labels (np.array): Indices of nodes forming tetrahedra
        
        Returns:
        ---
        centres (np.array): Array of radii and centres (shape: n_tetra x 4)
        """
        centres = np.zeros((len(cell_labels), 4))
        for i, n in enumerate(cell_labels):
            C, r = self.tetrahedron_circumscribed_sphere(nodes[n, :].T)
            centres[i] = [r, *C]
        return centres

    def load_global_files(self):
        """
        Load global mesh coordinate, tetrahedra index and neighbour files.
        """   
        self.global_mesh_node_coords = self.read_txt(self.global_mesh_node_coords, dtype=float) # Node coords
        self.global_mesh_tetra_indices = self.read_txt(self.global_mesh_tetra_indices, dtype=int, index_file=True) # Tetrahedra node indices
        self.global_mesh_tetra_neighbours = self.read_txt(self.global_mesh_tetra_neighbours, dtype=int, index_file=True) # Shared tetrahedra faces

        return

    def match_faces(self):
        """
        Map inner and outer surface face nodes to global mesh node indices.
        """
        # Load surface files
        self.out_surface_nodes = self.read_txt(os.path.join(self.surface_dir, "outer_surface_node_coords.txt"), dtype=float)
        self.out_surface_faces = self.read_txt(os.path.join(self.surface_dir, "outer_surface_face_indices.txt"), dtype=int, index_file=True)
        self.in_surface_nodes = self.read_txt(os.path.join(self.surface_dir, "inner_surface_node_coords.txt"), dtype=float)
        self.in_surface_faces = self.read_txt(os.path.join(self.surface_dir, "inner_surface_face_indices.txt"), dtype=int, index_file=True)
        
        inputs = {
            "outer" : [self.out_surface_nodes, self.out_surface_faces],
            "inner" : [self.in_surface_nodes, self.in_surface_faces]
        }
        
        for face in inputs:
            # Build KDTree for nearest neighbour matching
            tree = KDTree(self.global_mesh_node_coords)
            dists, matched = tree.query(inputs[face][0])

            # Check node matching against tolerances
            max_tol = 1e-3
            if np.any(dists > max_tol):
                self.helpers.plugin_log(f"WARNING: Some {face} surface nodes mapped > {max_tol}. Max dist = {dists.max():.4e}")

            # Remap face indices to global indices
            remapped = matched[inputs[face][1]]
    
            # Save output file
            output_file = os.path.join(self.interim_dir, f"{face}_surface_face_indices_mapped.txt")
            self.save_txt(output_file, remapped, index_file=True, add_index_col=False)
    
            # Check file produced
            if not os.path.isfile(output_file):
                self.helpers.errors(f"Mapping surface face file file not produced - {output_file}")

        return
        
    def write_mesh_info(self):
        """
        Write a txt file containing all mesh information.
        """
        # Determine indexing system
        self.idx_offset = 0 if self.parameters["zero_index_outputs"] else 1
        
        # Load files
        self.outer_faces = self.read_txt(os.path.join(self.interim_dir, f"outer_surface_face_indices_mapped.txt"), dtype=int, index_file=True)
        self.inner_faces = self.read_txt(os.path.join(self.interim_dir, f"inner_surface_face_indices_mapped.txt"), dtype=int, index_file=True)
        
        # Writing function
        def write_block(f, header, elems, prefix):
            f.write(f"{header}\n{len(elems)}\n")
            for e in elems:
                f.write(prefix.format(*e))

        # Write mesh file
        output_file = os.path.join(self.output_dir, "mesh_info.txt")
        with open(output_file, 'w') as f:
            f.write('$Node\n')
            f.write(f'{len(self.global_mesh_node_coords)}\n')
            nodes = np.column_stack([
                np.arange(self.idx_offset, len(self.global_mesh_node_coords) + self.idx_offset),
                self.global_mesh_node_coords
            ])
            for row in nodes:
                f.write(f"{int(row[0])} {row[1]:.6f} {row[2]:.6f} {row[3]:.6f}\n")

            # Save info to file   
            write_block(f, "$OuterFaceCell", self.outer_faces + self.idx_offset, "3 o {0} {1} {2}\n")
            write_block(f, "$InnerFaceCell", self.inner_faces + self.idx_offset, "3 i {0} {1} {2}\n")
            write_block(f, "$TetraCell", self.global_mesh_tetra_indices + self.idx_offset, "4 {0} {1} {2} {3}\n")

        # Check file produced
        if not os.path.isfile(output_file):
            self.helpers.errors(f"Mesh information file not produced - {output_file}")

        return

    def classify_tetrahedra(self):
        """
        Classify which region each tetrahedral element belongs
        to based on circumsphere and centres
        """
        # Define regions
        regions = {
            "cerebrum": [1, 2],
            "WM": [3, 4],
            "brainstem": [5, 6],
            "cerebellum": [7, 8],
            "cerebellumWM": [9, 10],
        }
    
        # Copmute centroids of global tetras
        tetra_centre = np.column_stack([
            self.global_mesh_node_coords[self.global_mesh_tetra_indices, :].mean(axis=1)
        ])

        # Loop over each region
        for region in regions:
            # Load region-specific coords and indices
            lh_node_coords   = self.read_txt(os.path.join(self.regional_mesh_dir, region, f"{region}_L_node_coords.txt"), dtype=float)
            lh_tetra_indices = self.read_txt(os.path.join(self.regional_mesh_dir, region, f"{region}_L_tetra_indices.txt"), dtype=int, index_file=True)
            rh_node_coords   = self.read_txt(os.path.join(self.regional_mesh_dir, region, f"{region}_R_node_coords.txt"), dtype=float)
            rh_tetra_indices = self.read_txt(os.path.join(self.regional_mesh_dir, region, f"{region}_R_tetra_indices.txt"), dtype=int, index_file=True)
    
            # Compute circumspheres (radius and centre)
            centres_lh = self.compute_tetra_centres(lh_node_coords, lh_tetra_indices)
            centres_rh = self.compute_tetra_centres(rh_node_coords, rh_tetra_indices)

            # Build KDTrees for nearest-neighbour search of centres
            tree_lh, tree_rh = KDTree(centres_lh[:, 1:]), KDTree(centres_rh[:, 1:])
            # Store maximum radius per hemisphere
            max_r_lh, max_r_rh = centres_lh[:, 0].max(), centres_rh[:, 0].max()
    
            # Initialise region labels
            labels = np.zeros((len(tetra_centre), 1), dtype=int)

            # Loop global tetrahedron centroids
            for i, pt in enumerate(tetra_centre):
                # Check left and right regional meshes
                for centres, tree, r_max, label in [
                    (centres_lh, tree_lh, max_r_lh, regions[region][0]),
                    (centres_rh, tree_rh, max_r_rh, regions[region][1]),
                ]:
                    # Skip if already labelled
                    if labels[i, 0]:
                        continue

                    # Find local cetnres within max radius of current point
                    idxs = tree.query_ball_point(pt, r=r_max)

                    # Check if point lies in any circumsphere
                    if idxs and np.any(np.sum((centres[idxs, 1:] - pt) ** 2, 1) < centres[idxs, 0] ** 2):
                        labels[i, 0] = label
    
            # Save labels (index + label)
            outpath = os.path.join(self.interim_dir, f"{region}_labels.txt")
            self.save_txt(outpath, labels)
    
        # Combine all regions, starting from cerebrum
        combined = self.read_txt(os.path.join(self.interim_dir, "cerebrum_labels.txt"), dtype=int, remove_idx_col=False)
        for region in regions:
            if region != "cerebrum":
                # Load new region label
                new = self.read_txt(os.path.join(self.interim_dir, f"{region}_labels.txt"), dtype=int)
                vals = new[:, 0].astype(int) # Pull vals from file
                idx = np.arange(len(vals)) # Create indices
                mask = (vals != 0) | (combined[idx, 1] == 0) # Update if new label is non-zero or existing is zero
                combined[idx[mask], 1] = vals[mask] # Apply new labels
    
        # Save combined labels
        output_file = os.path.join(self.interim_dir, "regional_labels.txt")
        self.save_txt(output_file, combined, add_index_col=False)
        self.labels_file = output_file

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
            "tensor": os.path.join(self.dwi_dir, "dwi_tensor.nii"),
            "L1": os.path.join(self.dwi_dir, "dwi_L1.nii"),
            "L2": os.path.join(self.dwi_dir, "dwi_L2.nii"),
            "L3": os.path.join(self.dwi_dir, "dwi_L3.nii"),
            "FA": os.path.join(self.dwi_dir, "dwi_FA.nii"),
            "MD": os.path.join(self.dwi_dir, "dwi_MD.nii"),
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

        # Read global mesh files
        self.load_global_files()

        # Mapping inner and outer surface nodes and faces to global indices
        self.helpers.plugin_log(f"Mapping inner and outer nodes and faces to global indices")
        self.match_faces()

        # Write mesh information file
        if self.parameters["write_mesh_info"]:
            self.helpers.plugin_log(f"Writing mesh information file")
            self.write_mesh_info()

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
        self.helpers.plugin_log(f"Core analysis completed at {self.helpers.now_time()}...")

        return