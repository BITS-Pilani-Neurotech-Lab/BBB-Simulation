# setup
bg_color black
hide everything

# load structure
load data/structures/complex_30A.pdb, complex

# receptor represented in gray
show cartoon, complex
color gray50, complex
set cartoon_transparency, 0.3, complex

# map and style fusion ligand
select fusion_construct, complex and (chain C or chain B or resi 1-145)
show cartoon, fusion_construct

color slate, fusion_construct and resi 1-71        # aga2-p scaffold
color yellow, fusion_construct and resi 73-109     # linker
color red, fusion_construct and resi 110-145       # rvg

# highlight close contacts
select clash_contacts, (complex and not fusion_construct) within 2.8 of fusion_construct
select clash_sidechains, clash_contacts or (fusion_construct within 2.8 of (complex and not fusion_construct))

show sticks, clash_sidechains
color hotpink, clash_sidechains
util.cnc clash_sidechains
set stick_radius, 0.25, clash_sidechains

# show bump distance markers
distance bump_dashes, complex and not fusion_construct, fusion_construct, 2.8, mode=2
color red, bump_dashes
set dash_width, 3.5
set dash_gap, 0.15

# quality
set cartoon_fancy_helices, 1
set antialias, 2

# save
draw 1920, 1080
png results/figures/3d_steric_clashes_docking_strain.png, dpi=300