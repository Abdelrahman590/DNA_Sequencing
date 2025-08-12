import sys
import os
from analyzer import DNAAnalyzer

def print_dict(d, indent=0):
    for key, value in d.items():
        if isinstance(value, dict):
            print(" " * indent + f"{key}:")
            print_dict(value, indent + 2)
        else:
            print(" " * indent + f"{key}: {value}")

def main():
    # Ensure proper encoding for Arabic text
    if sys.stdout.encoding != 'utf-8':
        sys.stdout.reconfigure(encoding='utf-8')

    # Test file path
    test_file = "../Data/chimpanzee.txt/chimpanzee.txt"
    abs_path = os.path.abspath(test_file)
    
    print(f"Testing file: {abs_path}")
    if not os.path.exists(test_file):
        print("Error: File not found!")
        return

    # Create analyzer instance
    analyzer = DNAAnalyzer()

    try:
        print("\nAnalyzing file...")
        results = analyzer.analyze_from_file(test_file)
        
        print("\nAnalysis Results:")
        if results:
            print_dict(results)
        else:
            print("No results returned")
            
    except Exception as e:
        print(f"\nError during analysis: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
