#!/usr/bin/env python 
# coding: utf-8

## Generates groups catalogs and stores them in `precomputed/` as json files.
## The keys are fingerprints of the groups.
## The values are the names of the groups.

load("main_functions.sage")

def main():
    ## input
    should_generate_tubical = True
    tubical_orderbound = 2000
    should_generate_toroidal = True
    toroidal_orderbound = 130
    should_generate_polyhedral_and_axial = True

    ## output directory
    output_dir = "precomputed/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    ## tubical
    if should_generate_tubical:
        catalog = {}
        generate_tubical(orderbound=tubical_orderbound, generate_full=False, names_only=True, verbose=False)
        assert len(set(flatten(list(catalog.values())))) == len(flatten(list(catalog.values())))
        filename = f"tubical_catalog_orderbound_{tubical_orderbound}.json"
        with open(output_dir + filename, "w") as f:
            json.dump(catalog, f, indent=4)
        print()

    ## toroidal
    if should_generate_toroidal:
        catalog = {}
        generate_chiral_toroidal(orderbound=toroidal_orderbound, generate_full=False, names_only=True, verbose=False)
        assert len(set(flatten(list(catalog.values())))) == len(flatten(list(catalog.values())))
        filename = f"toroidal_chiral_catalog_orderbound_{toroidal_orderbound}.json"
        with open(output_dir + filename, "w") as f:
            json.dump(catalog, f, indent=4)
        print()

        catalog = {}
        generate_achiral_toroidal(orderbound=toroidal_orderbound, generate_full=False, names_only=True, verbose=False)
        assert len(set(flatten(list(catalog.values())))) == len(flatten(list(catalog.values())))
        filename = f"toroidal_achiral_catalog_orderbound_{toroidal_orderbound}.json"
        with open(output_dir + filename, "w") as f:
            json.dump(catalog, f, indent=4)
        print()

    ## polyhedral and axial
    if should_generate_polyhedral_and_axial:
        catalog = {}
        generate_polyhedral_and_axial(names_only=True, verbose=False)
        filename = f"polyhedral_and_axial_catalog.json"
        with open(output_dir + filename, "w") as f:
            json.dump(catalog, f, indent=4)
        print()

if __name__ == "__main__":
    main()
