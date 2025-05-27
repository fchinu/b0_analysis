import csv
import sys

def read_numbers_from_csv(filename):
    """Read all numbers from a CSV file and return as a set of integers."""
    numbers = set()
    try:
        with open(filename, 'r', newline='', encoding='utf-8') as file:
            csv_reader = csv.reader(file)
            for row in csv_reader:
                for cell in row:
                    # Strip whitespace and skip empty cells
                    cell = cell.strip()
                    if cell:
                        try:
                            numbers.add(int(cell))
                        except ValueError:
                            print(f"Warning: '{cell}' is not a valid integer, skipping...")
        return numbers
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return set()
    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return set()

def find_common_numbers(file1, file2):
    """Find numbers that are present in both CSV files."""
    print(f"Reading numbers from {file1}...")
    numbers1 = read_numbers_from_csv(file1)
    print(f"Found {len(numbers1)} unique numbers in {file1}")
    
    print(f"Reading numbers from {file2}...")
    numbers2 = read_numbers_from_csv(file2)
    print(f"Found {len(numbers2)} unique numbers in {file2}")
    
    # Find intersection (common numbers)
    common = numbers1.intersection(numbers2)
    
    return sorted(common)

def main():
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <file1.csv> <file2.csv>")
        print("Example: python script.py numbers1.csv numbers2.csv")
        return
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    print("Finding common numbers between CSV files...")
    print("=" * 50)
    
    common_numbers = find_common_numbers(file1, file2)
    
    print("=" * 50)
    if common_numbers:
        print(f"Found {len(common_numbers)} common numbers:")
        print(common_numbers)
        
        output_file = "common_numbers.csv"
        with open(output_file, 'w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow(common_numbers)
        print(f"\nResults saved to: {output_file}")
    else:
        print("No common numbers found between the two files.")

if __name__ == "__main__":
    main()