import sys
from astropy.table import Table

def inspect(path):
    print(f"Inspecting: {path}")
    try:
        t = Table.read(path, format='fits')
        print(f"Rows: {len(t)}")
        print("Columns found:")
        print(t.colnames)
    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    inspect(sys.argv[1])
