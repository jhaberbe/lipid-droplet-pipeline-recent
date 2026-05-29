import argparse
import anndata

def main():
    parser = argparse.ArgumentParser(description="Single Cell Analysis Pipeline.")
    parser.add_argument("adata_path", type=str, help="Path to the .h5ad file")

    args = parser.parse_args()
    
    print(f"Loading AnnData from: {args.adata_path}")
    adata = anndata.read_h5ad(args.adata_path)
    
    print(f"AnnData loaded with shape: {adata.shape}")




if __name__ == "__main__":
    main()
