import numpy as np
from pymatgen.io.cif import CifParser
import warnings

def extract_cif_data_robust(cif_file_path):
    """
    Extract crystallographic data from CIF file with robust error handling
    """
    
    # Suppress warnings for cleaner output
    warnings.filterwarnings('ignore')
    
    try:
        # Parse the CIF file
        parser = CifParser(cif_file_path)
        structures = parser.parse_structures(primitive=False)
        structure = structures[0]
        
        # Initialize results dictionary
        results = {}
        
        # Basic Structure Information
        results['formula'] = str(structure.formula)
        results['reduced_formula'] = str(structure.reduced_formula)
        results['composition'] = dict(structure.composition.as_dict())
        results['num_sites'] = structure.num_sites
        results['volume'] = float(structure.volume)
        results['density'] = float(structure.density)
        
        # Lattice Parameters
        lattice = structure.lattice
        results['lattice'] = {
            'a': float(lattice.a),
            'b': float(lattice.b),
            'c': float(lattice.c),
            'alpha': float(lattice.alpha),
            'beta': float(lattice.beta),
            'gamma': float(lattice.gamma),
            'volume': float(lattice.volume)
        }
        
        # Crystal System and Space Group (with error handling)
        try:
            space_group_info = structure.get_space_group_info()
            if len(space_group_info) >= 2:
                results['crystal_system'] = str(space_group_info[0])
                results['space_group_number'] = int(space_group_info[1])
            else:
                results['crystal_system'] = "Unknown"
                results['space_group_number'] = "Unknown"
        except:
            results['crystal_system'] = "Unknown"
            results['space_group_number'] = "Unknown"
        
        # Atomic Information
        results['sites'] = []
        for i, site in enumerate(structure.sites):
            site_info = {
                'index': i,
                'species': str(site.specie),
                'fractional_coords': [float(x) for x in site.frac_coords],
                'cartesian_coords': [float(x) for x in site.coords],
                'occupancy': float(getattr(site, 'occupancy', 1.0))
            }
            results['sites'].append(site_info)
        
        # Simplified bond analysis - only calculate distances between atoms in unit cell
        results['bonds'] = []
        try:
            # Calculate distances between atoms within the unit cell only
            for i in range(len(structure.sites)):
                for j in range(i+1, len(structure.sites)):
                    site1 = structure.sites[i]
                    site2 = structure.sites[j]
                    
                    # Calculate distance
                    distance = structure.get_distance(i, j)
                    
                    # Only include reasonable bond distances (< 4.0 Å)
                    if distance < 4.0:
                        bond_info = {
                            'atom1_index': i,
                            'atom1_species': str(site1.specie),
                            'atom2_index': j,
                            'atom2_species': str(site2.specie),
                            'distance': float(distance)
                        }
                        results['bonds'].append(bond_info)
        except Exception as e:
            print(f"Note: Bond analysis skipped due to: {e}")
            results['bonds'] = []
        
        return results, structure
        
    except Exception as e:
        print(f"Error processing CIF file: {e}")
        return None, None

def print_cif_summary_safe(results):
    """
    Print a formatted summary with safe ASCII output
    """
    if not results:
        print("No data to display")
        return
        
    print("=" * 60)
    print("CRYSTALLOGRAPHIC DATA SUMMARY")
    print("=" * 60)
    
    print(f"Formula: {results['formula']}")
    print(f"Reduced Formula: {results['reduced_formula']}")
    print(f"Number of atoms: {results['num_sites']}")
    print(f"Volume: {results['volume']:.3f} A^3")
    print(f"Density: {results['density']:.3f} g/cm^3")
    
    print("\nLATTICE PARAMETERS:")
    print("-" * 30)
    lattice = results['lattice']
    print(f"a = {lattice['a']:.4f} A")
    print(f"b = {lattice['b']:.4f} A")
    print(f"c = {lattice['c']:.4f} A")
    print(f"alpha = {lattice['alpha']:.2f} degrees")
    print(f"beta = {lattice['beta']:.2f} degrees")
    print(f"gamma = {lattice['gamma']:.2f} degrees")
    
    print(f"\nCRYSTAL SYSTEM: {results['crystal_system']}")
    print(f"SPACE GROUP: {results['space_group_number']}")
    
    print("\nATOMIC POSITIONS (first 20 atoms):")
    print("-" * 60)
    print("Index | Species | Fractional Coordinates")
    print("-" * 60)
    for i, site in enumerate(results['sites']):  # Show first 20 atoms
        coords = site['fractional_coords']
        print(f"{site['index']:5d} | {site['species']:7s} | "
              f"({coords[0]:7.4f}, {coords[1]:7.4f}, {coords[2]:7.4f})")
    
    # Show composition breakdown
    print(f"\nCOMPOSITION:")
    print("-" * 30)
    for element, count in results['composition'].items():
        print(f"{element}: {count}")
    
    # Show some bond statistics if available
    if results['bonds']:
        print(f"\nBOND STATISTICS:")
        print("-" * 30)
        distances = [bond['distance'] for bond in results['bonds']]
        print(f"Total bonds (< 4.0 A): {len(distances)}")
        print(f"Average bond length: {np.mean(distances):.3f} A")
        print(f"Min bond length: {np.min(distances):.3f} A")
        print(f"Max bond length: {np.max(distances):.3f} A")
        
        print(f"\nSHORTEST BONDS (first 10):")
        print("-" * 50)
        sorted_bonds = sorted(results['bonds'], key=lambda x: x['distance'])
        for bond in sorted_bonds[:10]:
            print(f"{bond['atom1_species']}-{bond['atom2_species']}: {bond['distance']:.3f} A")

def analyze_composition(results):
    """
    Analyze the chemical composition
    """
    if not results:
        return
        
    composition = results['composition']
    total_atoms = sum(composition.values())
    
    print(f"\nDETAILED COMPOSITION ANALYSIS:")
    print("-" * 40)
    print(f"Total atoms: {total_atoms}")
    
    for element, count in sorted(composition.items()):
        percentage = (count / total_atoms) * 100
        print(f"{element}: {count} atoms ({percentage:.1f}%)")

def get_element_specific_bonds(results, element1, element2=None):
    """
    Get bonds involving specific elements
    """
    if not results or not results['bonds']:
        return []
    
    specific_bonds = []
    for bond in results['bonds']:
        if element2 is None:
            # Any bond involving element1
            if element1 in [bond['atom1_species'], bond['atom2_species']]:
                specific_bonds.append(bond)
        else:
            # Specific bond between element1 and element2
            if ((bond['atom1_species'] == element1 and bond['atom2_species'] == element2) or
                (bond['atom1_species'] == element2 and bond['atom2_species'] == element1)):
                specific_bonds.append(bond)
    
    return specific_bonds

# Example usage
if __name__ == "__main__":
    # Replace with your CIF file path
    cif_file = "mof.cif"
    
    # Extract data
    data, structure = extract_cif_data_robust(cif_file)
    
    if data:
        # Print summary
        print_cif_summary_safe(data)
        
        # Analyze composition
        analyze_composition(data)
        
        # Look at specific bonds (example: Zn-O bonds)
        if 'Zn' in data['composition'] and 'O' in data['composition']:
            zn_o_bonds = get_element_specific_bonds(data, 'Zn', 'O')
            if zn_o_bonds:
                print(f"\nZn-O BONDS:")
                print("-" * 30)
                for bond in zn_o_bonds:
                    print(f"Zn-O: {bond['distance']:.3f} A")    
    else:
        print("Failed to extract data from CIF file")

# Additional utility function for coordination analysis
def analyze_coordination_environment(results, central_atom='Zn', max_distance=3.0):
    """
    Analyze coordination environment around a specific atom type
    """
    if not results:
        return
    
    print(f"\nCOORDINATION ANALYSIS FOR {central_atom}:")
    print("-" * 50)
    
    central_sites = [site for site in results['sites'] if site['species'] == central_atom]
    
    for site in central_sites:
        coordinating_atoms = []
        
        # Find bonds involving this central atom
        for bond in results['bonds']:
            if (bond['atom1_index'] == site['index'] and bond['distance'] <= max_distance):
                coordinating_atoms.append((bond['atom2_species'], bond['distance']))
            elif (bond['atom2_index'] == site['index'] and bond['distance'] <= max_distance):
                coordinating_atoms.append((bond['atom1_species'], bond['distance']))
        
        print(f"{central_atom}{site['index']:2d}: ", end="")
        if coordinating_atoms:
            coord_str = ", ".join([f"{species}({dist:.2f}A)" for species, dist in coordinating_atoms])
            print(f"[{len(coordinating_atoms)}] {coord_str}")
        else:
            print("No coordination found")