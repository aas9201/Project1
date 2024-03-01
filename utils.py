from plasmid import Part, Module, Plasmid

import re


def sort_sites(plasmid_sequence):
    """Sorts the imported plasmid sequence and does some structual validity checks"""
    predefined_sites = {
        "OriT": "CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA",
        "T0": "CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG",
        "T1": "CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT",
        "PshAI": "GACNNNNGTC",
        "SwaI": "ATTTAAAT",
        "AscI": "GGCGCGCC",
        "FseI": "GGCCGGCC",
        "PacI": "TTAATTAA",
        "SpeI": "ACTAGT",
        "SanDI": "GGGTCCC",
    }

    reference_point = plasmid_sequence.find(
        "TTAATTAA"
    )  # This just makes PacI as the starting reference point of every plasmid sequence due to it being circular
    if reference_point == -1:
        raise ValueError("PacI RE site not found")
    ordered_sequence = (
        plasmid_sequence[reference_point:] + plasmid_sequence[:reference_point]
    )

    correct_order = [
        "PacI",
        "SpeI",
        "T0",
        "SanDI",
        "SwaI",
        "PshAI",
        "OriT",
        "FseI",
        "AscI",
        "T1",
    ]
    site_positions = {}
    first_occurance = []

    for (
        site_name,
        site_sequence,
    ) in (
        predefined_sites.items()
    ):  # Subs in the "NNN" in the PshAI as this RE site could be any combination
        pattern = re.sub("N", "[ACGT]", site_sequence)
        positions = []

        for i in range(len(ordered_sequence)):
            if re.match(
                pattern, ordered_sequence[i:]
            ):  # Checks if these sites are within the plasmid sequence and what positions they are in
                positions.append(i + 1)

        if positions:
            site_positions[site_name] = positions
            first_occurance.append(
                (site_name, positions[0])
            )  # Creates a new dictionary to store site and position

        first_occurance.sort(
            key=lambda x: x[1]
        )  # Sorts them by where they first appeared, as an invalid plasmid may contain more of one of these sites which is invalid

    sorted_sites = {site: site_positions[site] for site, _ in first_occurance}

    for i in sorted_sites:
        if len(sorted_sites[i]) > 1:
            raise ValueError(
                f"The plasmid contains more than one: {i}"
            )  # Checks if there is more than one each non-variable site which is invalid

    plasmid_order = list(sorted_sites.keys())
    if (
        plasmid_order != correct_order
    ):  # Checks to see if the plasmid has correct order/structure
        raise ValueError("Please make sure plasmid is of valid SEVA Format")

    return sorted_sites


def get_module(plasmid_sequence, name, ID, module_type, description=""):
    """Creates a new module object by slicing the plasmid sequence depending on which module_type you want"""
    sorted_sites = sort_sites(
        plasmid_sequence
    )  # Sorts the plasmid so that it can be sliced easily

    if (
        module_type == "cargo"
    ):  # Slices the plasmid between restriction enzyme site to get Cargo module
        cargo_sequence = plasmid_sequence[
            (sorted_sites["PacI"][0] + 7) : sorted_sites["SpeI"][0] - 1
        ]
        cargo_part = Part("cargo_part", "cargo", cargo_sequence, "cargo")
        cargo_module = Module(name, ID, "cargo", description, [cargo_part])

        return cargo_module

    if (
        module_type == "marker"
    ):  # Slices the plasmid between restriction enzyme site to get Marker module
        marker_sequence = plasmid_sequence[
            (sorted_sites["SwaI"][0] + 7) : sorted_sites["PshAI"][0] - 1
        ]
        marker_part = Part("module_part", "module_part", marker_sequence, "abr")

        marker_module = Module(name, ID, "marker", description, fetched=True)
        marker_module.SwaI = Part("SwaI", "swai", "ATTTAAAT", "re")
        pshai = plasmid_sequence[
            (sorted_sites["PshAI"][0] - 1) : sorted_sites["PshAI"][0] + 9
        ]  # As the PshAI contains NNN this needs to be created dpending on imported sequence
        marker_module.PshAI = Part("PshAI", "pshai", pshai, "re")
        marker_module.structure = [marker_module.SwaI, marker_part, marker_module.PshAI]

        return marker_module

    if (
        module_type == "origin"
    ):  # Slices the plasmid between restriction enzyme site to get Origin module
        origin_sequence = plasmid_sequence[
            (sorted_sites["FseI"][0] + 7) : sorted_sites["AscI"][0] - 1
        ]
        origin_part = Part("origin_part", "origin_part", origin_sequence, "ori")
        origin_module = Module(name, ID, "origin", description, [origin_part])

        return origin_module


def get_scars(
    plasmid_sequence,
):  # Slices plasmid between certain sites to check for any potential scars as each SEVA plasmid has certain scars, some don't
    sorted_sites = sort_sites(plasmid_sequence)
    scars = {}
    check_scars = {
        "Scar after T1": (plasmid_sequence[0 : sorted_sites["PacI"][0] - 1]),
        "Scar before T0": (
            plasmid_sequence[(sorted_sites["SpeI"][0] + 5) : sorted_sites["T0"][0] - 1]
        ),
        "Scar after T0": (
            plasmid_sequence[
                (sorted_sites["SanDI"][0] + 6) : sorted_sites["SwaI"][0] - 1
            ]
        ),
        "Scar before OriT": (
            plasmid_sequence[
                (sorted_sites["PshAI"][0] + 9) : sorted_sites["OriT"][0] - 1
            ]
        ),
        "Scar after OriT": (
            plasmid_sequence[
                (sorted_sites["OriT"][0] + 245) : sorted_sites["FseI"][0] - 1
            ]
        ),
        "Scar before T1": (
            plasmid_sequence[
                (sorted_sites["AscI"][0] + 7) : sorted_sites["SwaI"][0] - 1
            ]
        ),
    }

    for key, value in check_scars.items():
        if value:
            scars[key] = value
    print(scars)
    return scars


def linear_construction(plasmid1, plasmid2, desired_plasmid):
    """A function which models the construction of a desired plasmid from two input plasmids
    The two input plasmids are of Type(Plasmid) and desired plasmid should be a list of three modules
    which are from the two input plasmids"""

    #Input validation
    if not isinstance(desired_plasmid, list):
        raise TypeError("Desired Plasmid must be of Type: list")
    for modules in desired_plasmid:
        if not isinstance(modules, Module):
            raise TypeError("List should contain Module types only")

        # Creates a list containing the module types and compares it with predefined structure to see if its correct
    correct_val = ["cargo", "marker", "origin"]
    order = []
    for module in desired_plasmid:
        order.append(module.module_type)
    if correct_val != order:
        raise ValueError(
            "Please ensure the modules are in the correct order: [Cargo,Marker,Origin]"
        )

    #Creates list of module unique_ids from each plasmid for identification
    desired_list = [module.unique_id for module in desired_plasmid]

    plasmid1_list = [module.unique_id for module in plasmid1.components]
    
    plasmid2_list = [module.unique_id for module in plasmid2.components]
   
    #Set creation to model which plasmid is going to be parent/donor
    desired_set = set(desired_list)
    plasmid1_set = set(plasmid1_list)
    plasmid2_set = set(plasmid2_list)

    #As only inserting one module into the plasmid, the desired plasmid must contain either one or two modules
    contribution1 = len(desired_set.intersection(plasmid1_set))
    contribution2 = len(desired_set.intersection(plasmid2_set))
    if contribution1 in (0,3) or contribution2 in (0,3):
        raise ValueError("Linear Construction is not possible")


    #See if plasmid1 will be the parent plasmid and if so, what module will need to be inserted into it
    if contribution1 > contribution2:
        parent_plasmid = plasmid1_list
        if parent_plasmid[0] != desired_list[0]:
            print(
                f"Cargo module from {plasmid2.name} has been inserted into the parent plasmid: {plasmid1.name}"
            )
        if parent_plasmid[1] != desired_list[1]:
            print(
                f"Marker module from {plasmid2.name} has been inserted into the parent plasmid: {plasmid1.name}"
            )
        if parent_plasmid[2] != desired_list[2]:
            print(
                f"Origin module from {plasmid2.name} has been inserted into the parent plasmid: {plasmid1.name}"
            )

    #See if plasmid2 will be the parent plasmid and if so, what module will need to be inserted into it
    if contribution1 < contribution2:
        parent_plasmid = plasmid2_list
        if parent_plasmid[0] != desired_list[0]:
            print(
                f"Cargo module from {plasmid1.name} has been inserted into the parent plasmid: {plasmid2.name}"
            )
        if parent_plasmid[1] != desired_list[1]:
            print(
                f"Marker module from {plasmid1.name} has been inserted into the parent plasmid: {plasmid2.name}"
            )
        if parent_plasmid[2] != desired_list[2]:
            print(
                f"Origin module from {plasmid1.name} has been inserted into the parent plasmid: {plasmid2.name}"
            )