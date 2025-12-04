import argparse
import sys
import numpy as np
import pandas as pd
import h5py
import gzip
import shutil
import os
import tempfile

def main():
    parser = argparse.ArgumentParser(description="Extract sub-matrix from gcmap file (HDF5 direct).")
    parser.add_argument("--gcmap", required=True, help="Path to .gcmap or .gcmap.gz file")
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g. 3R)")
    parser.add_argument("--start", type=int, required=True, help="Start bp")
    parser.add_argument("--end", type=int, required=True, help="End bp")
    parser.add_argument("--out", required=True, help="Output CSV path")
    
    args = parser.parse_args()

    if not os.path.exists(args.gcmap):
        print(f"Error: File not found: {args.gcmap}")
        sys.exit(1)

    # Handle decompression if necessary
    h5_path = args.gcmap
    temp_file = None
    
    if args.gcmap.endswith(".gz"):
        print(f"Decompressing {args.gcmap} to temporary file...")
        fd, temp_file = tempfile.mkstemp(suffix=".gcmap")
        os.close(fd)
        
        try:
            with gzip.open(args.gcmap, 'rb') as f_in:
                with open(temp_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            h5_path = temp_file
        except Exception as e:
            print(f"Error decompressing file: {e}")
            if os.path.exists(temp_file):
                os.remove(temp_file)
            sys.exit(1)

    try:
        with h5py.File(h5_path, 'r') as f:
            # List available chromosomes (top level keys)
            available_chroms = list(f.keys())
            print(f"Available chromosomes: {available_chroms}")
            
            # Resolve chromosome name
            target_chrom = args.chrom
            if target_chrom not in available_chroms:
                candidates = [f"chr{args.chrom}", args.chrom.replace("chr", "")]
                found = False
                for c in candidates:
                    if c in available_chroms:
                        target_chrom = c
                        found = True
                        break
                if not found:
                    print(f"Error: Chromosome {args.chrom} not found in file.")
                    sys.exit(1)
            
            print(f"Using chromosome: {target_chrom}")
            
            # Find resolution (dataset)
            # We look for '5kb' or similar key under the chromosome group
            grp = f[target_chrom]
            
            # Try to find the 5kb dataset
            target_res_key = None
            # Prioritize '5kb'
            if '5kb' in grp:
                target_res_key = '5kb'
            else:
                # List keys and look for resolution
                print(f"Keys in {target_chrom}: {list(grp.keys())}")
                # Fallback logic could go here
                sys.exit("Error: Could not find '5kb' resolution dataset.")

            dset = grp[target_res_key]
            
            # Get attributes
            bin_size = 5000 # Default
            if 'binsize' in dset.attrs:
                bin_size = dset.attrs['binsize']
                
            print(f"Bin size: {bin_size} bp")
            
            # Assuming chr starts at 0 (standard for these files based on inspection)
            chr_start_bp = 0 
            
            # Coordinate conversion
            rel_start = args.start - chr_start_bp
            rel_end = args.end - chr_start_bp
            
            if rel_start < 0: rel_start = 0
            
            idx_start = int(np.floor(rel_start / bin_size))
            idx_end   = int(np.ceil(rel_end / bin_size))
            
            print(f"Matrix shape: {dset.shape}")
            max_bin = dset.shape[0]
            
            if idx_end > max_bin:
                print(f"Warning: adjusted end bin from {idx_end} to {max_bin}")
                idx_end = max_bin

            print(f"Extracting bins [{idx_start}:{idx_end}]")
            
            # Slicing
            sub_matrix = dset[idx_start:idx_end, idx_start:idx_end]
            
            # Create DataFrame
            row_coords = [chr_start_bp + (i * bin_size) for i in range(idx_start, idx_end)]
            df = pd.DataFrame(sub_matrix, index=row_coords, columns=row_coords)
            
            os.makedirs(os.path.dirname(args.out), exist_ok=True)
            df.to_csv(args.out)
            print(f"Saved extracted matrix to {args.out}")

    except Exception as e:
        print(f"Error processing HDF5 file: {e}")
        sys.exit(1)
        
    finally:
        # Cleanup temp file
        if temp_file and os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                print("Cleaned up temporary file.")
            except OSError:
                pass

if __name__ == "__main__":
    main()
