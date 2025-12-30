#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys

def find_file(name, search_root="."):
    """Walks the directory tree to find a file."""
    for root, dirs, files in os.walk(search_root):
        if name in files:
            return os.path.join(root, name)
    return None

def main():
    print("="*60)
    print("AUTO-BUILDER: RESOLVING PATH CONFLICTS")
    print("="*60)

    # 1. Find the artifacts
    print(" -> Hunting for files...")
    tex_file = find_file("axis_final_report.tex", ".")
    img1 = find_file("fig_axes.png", ".")
    img2 = find_file("fig_separation.png", ".")

    if not tex_file:
        print("[ERROR] Could not find 'axis_final_report.tex'. Are you in the project root?")
        sys.exit(1)
    
    if not img1 or not img2:
        print("[ERROR] Could not find images. Did you run 'generate_plots.py'?")
        print(f"   fig_axes.png found? {img1 is not None}")
        print(f"   fig_separation.png found? {img2 is not None}")
        sys.exit(1)

    # 2. Unify locations (Move images to where the .tex file is)
    tex_dir = os.path.dirname(tex_file)
    print(f" -> Working Directory: {tex_dir}")
    
    # Copy images to the tex directory if they aren't already there
    if os.path.dirname(img1) != tex_dir:
        shutil.copy(img1, tex_dir)
        print(f" -> Moved {os.path.basename(img1)} to {tex_dir}")
        
    if os.path.dirname(img2) != tex_dir:
        shutil.copy(img2, tex_dir)
        print(f" -> Moved {os.path.basename(img2)} to {tex_dir}")

    # 3. Compile
    print(" -> Running PDFLaTeX (Twice for citations)...")
    
    # Run TWICE to fix citations
    # We use cwd=tex_dir to run the command INSIDE that folder
    try:
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "axis_final_report.tex"], 
                       cwd=tex_dir, check=True)
        print("    (Pass 1 Complete)")
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "axis_final_report.tex"], 
                       cwd=tex_dir, check=True)
        print("    (Pass 2 Complete)")
    except subprocess.CalledProcessError:
        print("\n[ERROR] LaTeX Compilation Failed.")
        print("Check 'axis_final_report.log' in the target directory for details.")
        sys.exit(1)

    print("\n" + "="*60)
    print(f"SUCCESS! PDF is ready at: {os.path.join(tex_dir, 'axis_final_report.pdf')}")
    print("="*60)

if __name__ == "__main__":
    main()
