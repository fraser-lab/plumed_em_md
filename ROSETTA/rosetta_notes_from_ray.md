Dealing with Lys acetylation (adapted from https://www.rosettacommons.org/demos/latest/public/patched_residue_types_in_centroid-fullatom_protocols/README):

move lys_acetylated.txt to database/chemical/residue_type_sets/centroid/patches/lys_acetylated.txt

add the line:
patches/lys_acetylated.txt
into
database/chemical/residue_type_sets/centroid/patches.txt

Edit ALY to LYS in Residue name! and make sure chains are A not 1A
remove Mg and GTP

1) Weights. There are two sampling schemes in the xml script: 1) CarteisanSampler and 2) full-atom relax in dual space (torsion and cartesian). You can think CartSampler is to sample large backbone conformational changes, and Relax is to slightly perturb the structure to drive it into an energy minimum.

My rationale behind higher density weight at Cartsampler is to make sure large backbone conformational sampling is guided by density. However, since density is not that informative around the residues we perturb, in the Relax stage (full-atom), we let Rosetta's full-atom energy function take the role to slightly perturb the loop, and reward conformations which are more protein-like.

For the comparison problem, you can simply rescore the structures with the same density weight and same scoring function.

2) If I understood correctly, you were worried that the neighboring residues might need to move to accommodate the newly sampled conformation in the loop. This part should be taken care by the energy minimizations in both CartSampler and Relax. The energy minimization is carried out for the whole protein, not just the sampled regions. Nevertheless, you can certainly play around adding movers into the protocol to see how things turn out.

-nstruct
-ignore_unrecognized_res
<Add mover="loaddens"/> - switched to command line


idealize_jd2.linuxgccrelease -database /netapp/home/jaimefraser/database -in::path ./ -in::file::s AC_SYMM.pdb  -no_optH -out::path ./ -out::path::pdb ./ -chainbreaks  -overwrite

rosetta_scripts.linuxgccrelease -database /netapp/home/jaimefraser/database -in::file::s AC_SYMM_ideal_edit.pdb -edensity::mapfile AC_SYMM_MASKED_10_ORIGIN_0.mrc -parser::protocol new_multi_ray.xml   -edensity::mapreso 3.5 -default_max_cycles 200 -edensity::cryoem_scatterers -out::suffix _asymm -crystal_refine -beta
