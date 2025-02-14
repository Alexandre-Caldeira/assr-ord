import glob
from pathlib import Path

def combine_matlab_files(input_directory='.', output_filename='combined_output.txt'):
    """
    Combines multiple MATLAB (.m) files into a single text file.
    
    Args:
        input_directory (str): Directory containing MATLAB files (default: current directory)
        output_filename (str): Name of the combined output file (default: 'combined_output.txt')
    """
    # Get list of all .m files in the specified directory
    matlab_files = glob.glob(f"{input_directory}/*.m")
    
    # Create full path for output file
    output_path = Path(output_filename)
    
    try:
        # Open output file for writing
        with open(output_path, 'w', encoding='utf-8') as outfile:
            # Process each MATLAB file
            for file_path in sorted(matlab_files):
                # Get just the filename without path
                filename = Path(file_path).name
                
                # Write marker line with filename
                outfile.write(f"\n=== Begin {filename} ===\n\n")
                
                # Read and write the file contents
                with open(file_path, 'r', encoding='utf-8') as infile:
                    outfile.write(infile.read())
                    
                # Add separator between files
                outfile.write("\n=== End " + filename + " ===\n\n")
        
        print(f"Successfully combined {len(matlab_files)} files into {output_path}")
        
    except Exception as e:
        print(f"Error processing files: {str(e)}")
        
# Example usage
if __name__ == "__main__":
    # Combine files in current directory
    combine_matlab_files()
    
    # Or specify a different directory and output file
    # combine_matlab_files('./matlab_files', 'all_matlab_code.txt')
