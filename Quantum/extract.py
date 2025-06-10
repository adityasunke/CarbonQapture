import numpy as np
from pymatgen.io.cif import CifParser
import warnings

def extract_cif_data_robust(cif_file_path):
    """
    Extract crystallographic data from CIF file with robust error handling
    """
    
    # Suppress warnings for cleaner output - prevents UserWarning and FutureWarning messages
    # that commonly occur during CIF parsing due to format inconsistencies
    warnings.filterwarnings('ignore')
    
    try:
        # Initialize the CIF parser with the file path
        parser = CifParser(cif_file_path)
        
        # Parse structures using the newer method (not deprecated get_structures)
        # primitive=False maintains the original unit cell without reducing to primitive cell
        structures = parser.parse_structures(primitive=False)
        
        # Most CIF files contain one structure, so we take the first one
        structure = structures[0]
        
        # Initialize dictionary to store all extracted data
        results = {}
        
        # Extract basic structural information from the pymatgen Structure object
        results['formula'] = str(structure.formula)  # Full chemical formula
        results['reduced_formula'] = str(structure.reduced_formula)  # Simplified formula
        results['composition'] = dict(structure.composition.as_dict())  # Element counts
        results['num_sites'] = structure.num_sites  # Total number of atomic sites
        results['volume'] = float(structure.volume)  # Unit cell volume in Ų
        results['density'] = float(structure.density)  # Calculated density in g/cm³
        
        # Extract lattice parameters from the structure's lattice object
        lattice = structure.lattice
        results['lattice'] = {
            'a': float(lattice.a),        # Lattice parameter a in Angstroms
            'b': float(lattice.b),        # Lattice parameter b in Angstroms
            'c': float(lattice.c),        # Lattice parameter c in Angstroms
            'alpha': float(lattice.alpha), # Angle alpha in degrees
            'beta': float(lattice.beta),   # Angle beta in degrees
            'gamma': float(lattice.gamma), # Angle gamma in degrees
            'volume': float(lattice.volume) # Lattice volume
        }
        
        # Extract crystal system and space group information with error handling
        # This section is wrapped in try-except because space group detection can fail
        try:
            space_group_info = structure.get_space_group_info()
            # Check if we got the expected tuple with at least 2 elements
            if len(space_group_info) >= 2:
                results['crystal_system'] = str(space_group_info[0])  # Crystal system (e.g., "triclinic")
                results['space_group_number'] = int(space_group_info[1])  # Space group number
            else:
                # Fallback if space group info is incomplete
                results['crystal_system'] = "Unknown"
                results['space_group_number'] = "Unknown"
        except:
            # If space group analysis fails completely, set as unknown
            results['crystal_system'] = "Unknown"
            results['space_group_number'] = "Unknown"
        
        # Extract atomic site information for each atom in the structure
        results['sites'] = []
        for i, site in enumerate(structure.sites):
            # Create a dictionary for each atomic site with comprehensive information
            site_info = {
                'index': i,  # Sequential index for referencing
                'species': str(site.specie),  # Element symbol (e.g., "Zn", "O", "C")
                'fractional_coords': [float(x) for x in site.frac_coords],  # Fractional coordinates (0-1)
                'cartesian_coords': [float(x) for x in site.coords],  # Cartesian coordinates in Angstroms
                'occupancy': float(getattr(site, 'occupancy', 1.0))  # Site occupancy (usually 1.0)
            }
            results['sites'].append(site_info)
        
        # Perform simplified bond analysis - calculate distances between atoms in unit cell
        results['bonds'] = []
        try:
            # Double loop to calculate distances between all pairs of atoms
            for i in range(len(structure.sites)):
                for j in range(i+1, len(structure.sites)):  # j > i to avoid duplicates
                    site1 = structure.sites[i]
                    site2 = structure.sites[j]
                    
                    # Calculate the distance between atoms i and j
                    # This accounts for periodic boundary conditions
                    distance = structure.get_distance(i, j)
                    
                    # Only store bonds with reasonable distances (< 4.0 Å)
                    # This filters out very long non-bonding interactions
                    if distance < 4.0:
                        bond_info = {
                            'atom1_index': i,  # Index of first atom
                            'atom1_species': str(site1.specie),  # Element of first atom
                            'atom2_index': j,  # Index of second atom
                            'atom2_species': str(site2.specie),  # Element of second atom
                            'distance': float(distance)  # Distance in Angstroms
                        }
                        results['bonds'].append(bond_info)
        except Exception as e:
            # If bond analysis fails, continue without bond information
            print(f"Note: Bond analysis skipped due to: {e}")
            results['bonds'] = []
        
        # Return both the processed data dictionary and the original structure object
        return results, structure
        
    except Exception as e:
        # If any part of the parsing fails, return None values
        print(f"Error processing CIF file: {e}")
        return None, None

def print_cif_summary_safe(results):
    """
    Print a formatted summary with safe ASCII output
    """
    # Check if we have valid results to display
    if not results:
        print("No data to display")
        return
        
    # Print header section with basic structural information
    print("=" * 60)
    print("CRYSTALLOGRAPHIC DATA SUMMARY")
    print("=" * 60)
    
    # Display fundamental composition and structural data
    print(f"Formula: {results['formula']}")
    print(f"Reduced Formula: {results['reduced_formula']}")
    print(f"Number of atoms: {results['num_sites']}")
    print(f"Volume: {results['volume']:.3f} A^3")  # Using A^3 instead of Ų for ASCII safety
    print(f"Density: {results['density']:.3f} g/cm^3")
    
    # Display lattice parameters section
    print("\nLATTICE PARAMETERS:")
    print("-" * 30)
    lattice = results['lattice']
    print(f"a = {lattice['a']:.4f} A")      # Unit cell dimension a
    print(f"b = {lattice['b']:.4f} A")      # Unit cell dimension b
    print(f"c = {lattice['c']:.4f} A")      # Unit cell dimension c
    print(f"alpha = {lattice['alpha']:.2f} degrees")  # Angle between b and c
    print(f"beta = {lattice['beta']:.2f} degrees")    # Angle between a and c
    print(f"gamma = {lattice['gamma']:.2f} degrees")  # Angle between a and b
    
    # Display symmetry information
    print(f"\nCRYSTAL SYSTEM: {results['crystal_system']}")
    print(f"SPACE GROUP: {results['space_group_number']}")
    
    # Display atomic positions - showing all atoms instead of limiting to 20
    print("\nATOMIC POSITIONS (first 20 atoms):")
    print("-" * 60)
    print("Index | Species | Fractional Coordinates")
    print("-" * 60)
    # Loop through all atomic sites and display their information
    for i, site in enumerate(results['sites']):  # Show first 20 atoms
        coords = site['fractional_coords']
        print(f"{site['index']:5d} | {site['species']:7s} | "
              f"({coords[0]:7.4f}, {coords[1]:7.4f}, {coords[2]:7.4f})")
    
    # Display composition breakdown showing element counts
    print(f"\nCOMPOSITION:")
    print("-" * 30)
    for element, count in results['composition'].items():
        print(f"{element}: {count}")
    
    # Display bond statistics if bond analysis was successful
    if results['bonds']:
        print(f"\nBOND STATISTICS:")
        print("-" * 30)
        # Extract all distances for statistical analysis
        distances = [bond['distance'] for bond in results['bonds']]
        print(f"Total bonds (< 4.0 A): {len(distances)}")
        print(f"Average bond length: {np.mean(distances):.3f} A")
        print(f"Min bond length: {np.min(distances):.3f} A")
        print(f"Max bond length: {np.max(distances):.3f} A")
        
        # Show the shortest bonds (most likely to be actual chemical bonds)
        print(f"\nSHORTEST BONDS (first 10):")
        print("-" * 50)
        sorted_bonds = sorted(results['bonds'], key=lambda x: x['distance'])
        for bond in sorted_bonds[:10]:
            print(f"{bond['atom1_species']}-{bond['atom2_species']}: {bond['distance']:.3f} A")

def analyze_composition(results):
    """
    Analyze the chemical composition
    """
    # Ensure we have valid results before proceeding
    if not results:
        return
        
    composition = results['composition']
    total_atoms = sum(composition.values())  # Calculate total number of atoms
    
    print(f"\nDETAILED COMPOSITION ANALYSIS:")
    print("-" * 40)
    print(f"Total atoms: {total_atoms}")
    
    # Display each element with count and percentage
    for element, count in sorted(composition.items()):
        percentage = (count / total_atoms) * 100  # Calculate percentage composition
        print(f"{element}: {count} atoms ({percentage:.1f}%)")

def get_element_specific_bonds(results, element1, element2=None):
    """
    Get bonds involving specific elements
    """
    # Return empty list if no results or no bonds available
    if not results or not results['bonds']:
        return []
    
    specific_bonds = []
    # Search through all bonds for those involving the specified elements
    for bond in results['bonds']:
        if element2 is None:
            # If only one element specified, find any bond involving that element
            if element1 in [bond['atom1_species'], bond['atom2_species']]:
                specific_bonds.append(bond)
        else:
            # If two elements specified, find bonds between those specific elements
            # Check both orientations since bonds are stored as atom1-atom2 pairs
            if ((bond['atom1_species'] == element1 and bond['atom2_species'] == element2) or
                (bond['atom1_species'] == element2 and bond['atom2_species'] == element1)):
                specific_bonds.append(bond)
    
    return specific_bonds

# Main execution block - runs when script is executed directly
if __name__ == "__main__":
    # Replace with your CIF file path
    cif_file = "mof.cif"
    
    # Extract data using the robust extraction function
    data, structure = extract_cif_data_robust(cif_file)
    
    # Process and display results if extraction was successful
    if data:
        # Print comprehensive summary of the structure
        print_cif_summary_safe(data)
        
        # Analyze and display composition percentages
        analyze_composition(data)
        
        # Look at specific bonds (example: Zn-O bonds) if both elements are present
        if 'Zn' in data['composition'] and 'O' in data['composition']:
            zn_o_bonds = get_element_specific_bonds(data, 'Zn', 'O')
            if zn_o_bonds:
                print(f"\nZn-O BONDS:")
                print("-" * 30)
                # Display all Zn-O bond distances
                for bond in zn_o_bonds:
                    print(f"Zn-O: {bond['distance']:.3f} A")    
    else:
        print("Failed to extract data from CIF file")

# Additional utility function for coordination analysis
def analyze_coordination_environment(results, central_atom='Zn', max_distance=3.0):
    """
    Analyze coordination environment around a specific atom type
    """
    # Ensure we have valid results before proceeding
    if not results:
        return
    
    print(f"\nCOORDINATION ANALYSIS FOR {central_atom}:")
    print("-" * 50)
    
    # Find all sites containing the central atom of interest
    central_sites = [site for site in results['sites'] if site['species'] == central_atom]
    
    # Analyze coordination environment for each central atom
    for site in central_sites:
        coordinating_atoms = []
        
        # Find all bonds involving this central atom within the specified distance
        for bond in results['bonds']:
            # Check if this central atom is involved in the bond and within distance limit
            if (bond['atom1_index'] == site['index'] and bond['distance'] <= max_distance):
                coordinating_atoms.append((bond['atom2_species'], bond['distance']))
            elif (bond['atom2_index'] == site['index'] and bond['distance'] <= max_distance):
                coordinating_atoms.append((bond['atom1_species'], bond['distance']))
        
        # Display coordination information for this central atom
        print(f"{central_atom}{site['index']:2d}: ", end="")
        if coordinating_atoms:
            # Create formatted string showing coordination number and coordinating atoms
            coord_str = ", ".join([f"{species}({dist:.2f}A)" for species, dist in coordinating_atoms])
            print(f"[{len(coordinating_atoms)}] {coord_str}")
        else:
            print("No coordination found")