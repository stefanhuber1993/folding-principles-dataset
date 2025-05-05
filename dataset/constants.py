# DSSP uses many different codes to represent secondary structure.
# This mapping simplifies them to a few categories H (helix), E (strand), C (Coil).
simplified_ss_map = {
    "H": "H", "I": "H",
    "E": "E", "B": "E",
    "G": "C", "T": "C", "S": "C", "-": "C"
}
